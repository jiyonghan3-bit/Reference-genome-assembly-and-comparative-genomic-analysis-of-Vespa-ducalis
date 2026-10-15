#!/usr/bin/env python3
"""Alignment and isoform integrity filters applied to every single-copy orthogroup
before the branch and branch-site tests.

Usage     python 07_filter_alignment_isoform.py [threads]
Requires  MAFFT (tested with 7.525), Python 3.8+
"""
import os
import sys
import csv
import subprocess
import tempfile
import statistics as st
from collections import Counter
from multiprocessing import Pool

# paths
MAFFT = "/PATH/TO/YOUR/TOOLS/mafft"
ALN_DIR = "/PATH/TO/YOUR/RESULTS/alignment/prank"              # from 05, untrimmed <orthogroup>.best.fas
USED_ISOFORMS = "/PATH/TO/YOUR/RESULTS/palo/palo_selected_isoforms.tsv"   # from 03
PROT2GENE = "/PATH/TO/YOUR/RESULTS/palo/tables/{species}.prot2gene.tsv"   # from 03
PROTEOME = "/PATH/TO/YOUR/RESULTS/proteomes/{species}.faa"     # EVERY isoform with its original
                                                               # id, not the longest or the PALO
                                                               # pick: the rescue checks need them
CODEML_SUMMARY = "/PATH/TO/YOUR/RESULTS/codeml_results.header.tsv"        # from 06, with column names
CODEML_DIR = "/PATH/TO/YOUR/RESULTS/codeml"                    # <orthogroup>/branch/mlc
OUT_DIR = "/PATH/TO/YOUR/RESULTS/filter"

COL_HOG = "HOG"
COL_DS = "fg_dS_M0"
COL_LEN = "ncodon"
COL_FGBR = "fg_branch"

SPECIES = ["Vduc", "Vman", "Vcra", "Vvel", "Vvul", "Vpen"]
FG = "Vduc"
BACKGROUND = [s for s in SPECIES if s != FG]

FG_DS_MIN = 0.0005   # dS, one-ratio model
BR_DS_MIN = 0.0      # dS, two-ratio model
LEN_TOL = 0.05       # (a)
BG_TOL = 1.05        # (a)
RESCUE_TOL = 0.30    # (c), (d)

CODON = {}
for _i, _b1 in enumerate("TCAG"):
    for _b2 in "TCAG":
        for _b3 in "TCAG":
            CODON[_b1 + _b2 + _b3] = (
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
            )[_i * 16 + "TCAG".index(_b2) * 4 + "TCAG".index(_b3)]

AMBIGUOUS = set("X*")


def read_fasta(path):
    seqs, key = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            key = line[1:].split()[0]
            seqs[key] = []
        elif key is not None:
            seqs[key].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def translate(seq):
    if len(seq) % 3:
        raise ValueError(f"alignment length {len(seq)} is not a multiple of three")
    return "".join(
        "-" if seq[i:i + 3] == "---" else CODON.get(seq[i:i + 3].upper(), "X")
        for i in range(0, len(seq), 3)
    )


def residue_index(seq):
    col2res, res2col = [], []
    n = 0
    for ch in seq:
        if ch == "-":
            col2res.append(None)
        else:
            col2res.append(n)
            res2col.append(len(col2res) - 1)
            n += 1
    return col2res, res2col


def differences(fg_seq, bg_seqs):
    # columns where all background species share one residue and the foreground differs
    out = []
    for i, residue in enumerate(fg_seq):
        if residue == "-" or residue in AMBIGUOUS:
            continue
        col = [s[i] for s in bg_seqs]
        if any(c == "-" or c in AMBIGUOUS for c in col):
            continue
        if len(set(col)) == 1 and col[0] != residue:
            out.append((i, residue, col[0]))
    return out


def mafft(seqs):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")
        tmp = fh.name
    try:
        proc = subprocess.run([MAFFT, "--auto", "--quiet", tmp],
                              capture_output=True, text=True)
    finally:
        os.unlink(tmp)
    if proc.returncode != 0 or not proc.stdout.strip():
        raise RuntimeError(f"MAFFT failed (exit {proc.returncode}): {proc.stderr[:200]}")
    aln, key = {}, None
    for line in proc.stdout.split("\n"):
        if line.startswith(">"):
            key = line[1:].split()[0]
            aln[key] = []
        elif line.strip() and key:
            aln[key].append(line.strip())
    return {k: "".join(v) for k, v in aln.items()}


def branch_foreground_dS(hog, fg_key):
    # foreground dS and S*dS from the two-ratio model
    path = f"{CODEML_DIR}/{hog}/branch/mlc"
    if not fg_key or not os.path.exists(path):
        return None, None
    lines = open(path).read().split("\n")
    head = [i for i, l in enumerate(lines) if l.strip().startswith("branch") and "dN/dS" in l]
    if not head:
        return None, None
    for line in lines[head[0] + 1:]:
        f = line.split()
        if len(f) >= 9 and f[0] == fg_key:
            return float(f[6]), float(f[8])
    return None, None


PROT, P2G, G2P, USED = {}, {}, {}, {}


def load():
    for sp in SPECIES:
        PROT[sp] = {k: v.rstrip("*")
                    for k, v in read_fasta(PROTEOME.format(species=sp)).items()}
        P2G[sp], G2P[sp] = {}, {}
        for line in open(PROT2GENE.format(species=sp)):
            prot, gene = line.rstrip("\n").split("\t")[:2]
            P2G[sp][prot] = gene
            G2P[sp].setdefault(gene, []).append(prot)
    for row in csv.DictReader(open(USED_ISOFORMS), delimiter="\t"):
        USED[row[COL_HOG]] = {sp: row[sp] for sp in SPECIES}


def rescue(hog, direction, diffs):
    # realign the aligned isoforms with every other isoform of the same genes
    # direction "own"     another isoform of the foreground gene carries the background residue
    #           "reverse" another isoform of a background gene carries the foreground residue
    used = USED.get(hog)
    if not used:
        return None, "no_isoform_table"
    if any(used[sp] not in PROT[sp] for sp in SPECIES):
        return None, "missing_protein"
    if any(used[sp] not in P2G[sp] for sp in SPECIES):
        return None, "missing_gene_mapping"
    if not diffs:
        return 0.0, "no_difference"

    seqs = {sp: PROT[sp][used[sp]] for sp in SPECIES}
    alt = {}
    for sp in ([FG] if direction == "own" else BACKGROUND):
        gene = P2G[sp][used[sp]]
        for prot in G2P[sp].get(gene, []):
            if prot != used[sp] and prot in PROT[sp]:
                alt[f"{sp}#ALT#{prot}"] = PROT[sp][prot]
    if not alt:
        return 0.0, "single_isoform"

    try:
        aln = mafft({**seqs, **alt})
    except RuntimeError:
        return None, "mafft_failed"
    if FG not in aln or any(s not in aln for s in BACKGROUND):
        return None, "mafft_incomplete"

    # locate the same foreground residues in the realignment
    _, res2col = residue_index(aln[FG])
    rescued = unmapped = 0
    for res_idx, fg_residue, bg_residue in diffs:
        if res_idx >= len(res2col):
            unmapped += 1
            continue
        col = res2col[res_idx]
        if direction == "own":
            hit = any(aln[k][col] == bg_residue for k in alt if k.startswith(f"{FG}#"))
        else:
            hit = any(aln[k][col] == fg_residue for k in alt if not k.startswith(f"{FG}#"))
        rescued += int(hit)
    n = len(diffs) - unmapped
    if n <= 0:
        return None, "unmapped"
    note = f"{rescued}/{n}" + (f" unmapped={unmapped}" if unmapped else "")
    return rescued / n, note


def one(hog):
    path = f"{ALN_DIR}/{hog}.best.fas"
    if not os.path.exists(path):
        return dict(HOG=hog, status="no_alignment")
    nt = read_fasta(path)
    if any(s not in nt for s in SPECIES):
        return dict(HOG=hog, status="missing_species")
    if len({len(nt[s]) for s in SPECIES}) != 1:
        return dict(HOG=hog, status="ragged_alignment")
    try:
        aa = {s: translate(nt[s]) for s in SPECIES}
    except ValueError:
        return dict(HOG=hog, status="length_not_multiple_of_three")

    fg = aa[FG]
    aln_len = len(fg)

    # (a) length ratio
    lens = {s: len(aa[s].replace("-", "")) for s in SPECIES}
    bg = [lens[s] for s in BACKGROUND]
    median = st.median(bg)
    ratio = lens[FG] / median if median else 0
    bg_maxmin = max(bg) / min(bg) if min(bg) else float("inf")

    diffs_col = differences(fg, [aa[s] for s in BACKGROUND])
    n = len(diffs_col)

    # (c), (d) rescue
    col2res, _ = residue_index(fg)
    diffs_res = [(col2res[i], t, c) for i, t, c in diffs_col]
    c_rescue, c_note = rescue(hog, "own", diffs_res)
    d_rescue, d_note = rescue(hog, "reverse", diffs_res)

    return dict(
        HOG=hog, status="OK", aln_aa=aln_len,
        len_fg=lens[FG], len_bg_median=median, len_ratio=ratio, bg_maxmin=bg_maxmin,
        n_diff=n,
        self_rescue=c_rescue, self_rescue_n=c_note,
        reverse_rescue=d_rescue, reverse_rescue_n=d_note,
        pass_a=int(abs(ratio - 1) <= LEN_TOL and bg_maxmin <= BG_TOL),
        pass_c=int(c_rescue is not None and c_rescue < RESCUE_TOL),
        pass_d=int(d_rescue is not None and d_rescue < RESCUE_TOL),
    )


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    os.makedirs(OUT_DIR, exist_ok=True)

    codeml = {r[COL_HOG]: r for r in csv.DictReader(open(CODEML_SUMMARY), delimiter="\t")}
    branch_dS = {h: branch_foreground_dS(h, codeml[h].get(COL_FGBR, "")) for h in codeml}
    hogs = sorted(codeml)
    print(f"{len(hogs)} orthogroups, {threads} threads", flush=True)

    # initializer, not inherited globals: fork is not the default start method everywhere
    with Pool(threads, initializer=load) as pool:
        results = pool.map(one, hogs, chunksize=20)

    cols = ["HOG", "status", "codon_len", "aln_aa", "fg_dS", "fg_dS_branch",
            "fg_syn_subs", "pass_dS",
            "len_fg", "len_bg_median", "len_ratio", "bg_maxmin", "pass_a",
            "n_diff",
            "self_rescue", "self_rescue_n", "pass_c",
            "reverse_rescue", "reverse_rescue_n", "pass_d", "pass_all"]
    counts = Counter()
    states = Counter()

    with open(f"{OUT_DIR}/filter_results.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in results:
            states[row["status"]] += 1
            src = codeml[row["HOG"]]
            ds = float(src[COL_DS]) if src[COL_DS] else 0.0
            row["fg_dS"] = ds
            row["codon_len"] = src[COL_LEN]
            bds, bsubs = branch_dS.get(row["HOG"], (None, None))
            row["fg_dS_branch"] = "" if bds is None else bds
            row["fg_syn_subs"] = "" if bsubs is None else bsubs
            row["pass_dS"] = int(ds > FG_DS_MIN and bds is not None and bds > BR_DS_MIN)
            if row["status"] == "OK":
                gates = ("pass_dS", "pass_a", "pass_c", "pass_d")
                row["pass_all"] = int(all(row[k] for k in gates))
                for k in gates + ("pass_all",):
                    counts[k] += row[k]
            else:
                row["pass_all"] = 0
            counts["total"] += 1
            writer.writerow(row)

    summary = [
        "alignment and isoform integrity filters",
        f"  orthogroups examined                       {counts['total']}",
        f"  dS  foreground dS (M0 > {FG_DS_MIN}, branch > {BR_DS_MIN})   {counts['pass_dS']}",
        f"  (a) length ratio                           {counts['pass_a']}",
        f"  (c) self rescue                            {counts['pass_c']}",
        f"  (d) reverse rescue                         {counts['pass_d']}",
        f"  passed all gates (multiple-testing family) {counts['pass_all']}",
        "",
        "  status counts: " + ", ".join(f"{k}={v}" for k, v in sorted(states.items())),
    ]
    open(f"{OUT_DIR}/summary.txt", "w").write("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
