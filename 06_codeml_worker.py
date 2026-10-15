#!/usr/bin/env python
"""06_codeml_worker.py <orthogroup> — runs the four codeml models for one orthogroup
and prints a one-line summary. Called by 06_codeml.sh, one call per orthogroup.

Paths come from the environment: ALNDIR, TREE, OUTDIR, CODEML
Requires  PAML (tested with 4.10.10), Biopython
"""
import os, re, sys, time, shutil
from Bio.Phylo.PAML import codeml

HOG = sys.argv[1]
ALNDIR = os.environ["ALNDIR"]
TREE = os.environ["TREE"]
OUTDIR = os.environ["OUTDIR"]
CODEML = os.environ["CODEML"]

ALN = f"{ALNDIR}/{HOG}.fa"
D = f"{OUTDIR}/{HOG}"

COMMON = dict(seqtype=1, CodonFreq=2, clock=0, aaDist=0, icode=0, verbose=0,
              fix_kappa=0, kappa=2, fix_alpha=1, alpha=0, ncatG=1, getSE=0,
              RateAncestor=0, Small_Diff=.5e-6, cleandata=0, ndata=1)

# name -> (NSsites key, start values to try); the trial with the highest lnL is kept
RUNS = [
    ("M0",      0, [dict(model=0, NSsites=[0], fix_omega=0, omega=0.4),
                    dict(model=0, NSsites=[0], fix_omega=0, omega=1.5)]),
    ("branch",  0, [dict(model=2, NSsites=[0], fix_omega=0, omega=0.4),
                    dict(model=2, NSsites=[0], fix_omega=0, omega=1.5)]),
    ("bs_alt",  2, [dict(model=2, NSsites=[2], fix_omega=0, omega=0.5),
                    dict(model=2, NSsites=[2], fix_omega=0, omega=1.5),
                    dict(model=2, NSsites=[2], fix_omega=0, omega=0.5, kappa=4),
                    dict(model=2, NSsites=[2], fix_omega=0, omega=3.0)]),
    ("bs_null", 2, [dict(model=2, NSsites=[2], fix_omega=1, omega=1.0),
                    dict(model=2, NSsites=[2], fix_omega=1, omega=1.0, kappa=1),
                    dict(model=2, NSsites=[2], fix_omega=1, omega=1.0, kappa=4)]),
]
KEEP = {"mlc"}


def run_model(name, key, trials):
    best, best_dir, tried = None, f"{D}/{name}", []
    for i, opts in enumerate(trials):
        wd = best_dir if len(trials) == 1 else f"{D}/{name}.try{i}"
        os.makedirs(wd, exist_ok=True)
        cml = codeml.Codeml(alignment=ALN, tree=TREE, out_file=f"{wd}/mlc", working_dir=wd)
        cml.set_options(**{**COMMON, "method": 0, **opts})
        try:
            r = cml.run(command=CODEML, verbose=False, parse=True)
            ns = r["NSsites"][key]
            lnl = float(ns["lnL"])
        except Exception as e:
            sys.stderr.write(f"{HOG}\t{name}.try{i}\t{type(e).__name__}: {e}\n")
            tried.append("NA")
            continue
        tried.append(f"{lnl:.4f}")
        if best is None or lnl > best[0]:
            best = (lnl, ns, wd)
    if best is None:
        return None
    lnl, ns, wd = best
    if wd != best_dir:
        shutil.rmtree(best_dir, ignore_errors=True)
        shutil.move(wd, best_dir)
    for i in range(len(trials)):
        shutil.rmtree(f"{D}/{name}.try{i}", ignore_errors=True)
    for f in os.listdir(best_dir):
        if f not in KEEP:
            p = f"{best_dir}/{f}"
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
    return lnl, ns, ",".join(tried)


def fg_branch(ns, idx):
    # codeml numbers tips by their order in the alignment, so the branch key is read, not assumed
    for k, v in ns["parameters"].get("branches", {}).items():
        if k.endswith(f"..{idx}"):
            return k, v
    return "", {}


def parse_beb(mlc):
    try:
        txt = open(mlc).read()
    except OSError:
        return []
    m = re.search(r"Bayes Empirical Bayes \(BEB\).*?Prob\(w>1\):\s*\n(.*?)\n\s*\n", txt, re.S)
    if not m:
        return []
    out = []
    for line in m.group(1).split("\n"):
        mm = re.match(r"\s*(\d+)\s+(\S)\s+([\d.]+)", line)
        if mm:
            out.append((int(mm.group(1)), mm.group(2), float(mm.group(3))))
    return out


def main():
    t0 = time.time()
    os.makedirs(D, exist_ok=True)

    names = [l[1:].split()[0] for l in open(ALN) if l.startswith(">")]
    ncodon = sum(len(l.strip()) for l in open(ALN).read().split(">")[1].split("\n")[1:]) // 3
    if "Vduc" not in names:
        print(f"{HOG}\tFAIL\tforeground species missing", file=sys.stderr)
        return 1
    idx = names.index("Vduc") + 1

    res = {}
    for name, key, trials in RUNS:
        r = run_model(name, key, trials)
        if r is None:
            print(f"{HOG}\tFAIL\t{name}", file=sys.stderr)
            return 1
        res[name] = r

    lnl_m0, ns_m0, tr_m0 = res["M0"]
    fgk, fgb = fg_branch(ns_m0, idx)
    lnl_br, ns_br, tr_br = res["branch"]
    w = ns_br["parameters"]["omega"]
    w_bg, w_fg = (w[0], w[1]) if isinstance(w, (list, tuple)) else (w, w)
    lnl_ba, ns_ba, tr_ba = res["bs_alt"]
    lnl_bn, _, tr_bn = res["bs_null"]
    sc = ns_ba["parameters"].get("site classes", {})
    p2 = sum(sc[c]["proportion"] for c in (2, 3) if c in sc)
    w2 = sc.get(2, {}).get("branch types", {}).get("foreground", float("nan"))

    beb = parse_beb(f"{D}/bs_alt/mlc")
    with open(f"{D}/beb.tsv", "w") as f:
        for pos, aa, pr in beb:
            f.write(f"{HOG}\t{pos}\t{aa}\t{pr:.4f}\n")
    beb95 = [b for b in beb if b[2] > 0.95]

    print("\t".join(str(x) for x in [
        HOG, ncodon, fgk,
        f"{lnl_m0:.6f}", ns_m0["parameters"]["omega"],
        fgb.get("dS", ""), fgb.get("dN", ""), fgb.get("t", ""),
        f"{lnl_br:.6f}", w_bg, w_fg,
        f"{lnl_ba:.6f}", f"{lnl_bn:.6f}", f"{p2:.6f}", w2,
        len(beb), len(beb95),
        f"{max([b[2] for b in beb], default=0):.4f}",
        tr_m0, tr_br, tr_ba, tr_bn,
        f"{time.time() - t0:.1f}", "OK"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
