#!/usr/bin/env python3
"""Likelihood ratio tests and Benjamini-Hochberg correction.

Orthogroups are restricted to those that passed 07_filter_alignment_isoform.py, and the
correction is computed once over that set, separately for the branch and the
branch-site test.

Requires  scipy, statsmodels
"""
import os, csv, math
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

CODEML_SUMMARY = "/PATH/TO/YOUR/RESULTS/codeml_results.header.tsv"   # the file with column names
FILTER_RESULTS = "/PATH/TO/YOUR/RESULTS/filter/filter_results.tsv"
OUT = "/PATH/TO/YOUR/RESULTS/tests"
FDR_Q = 0.05
BEB_P = 0.95


def f(x):
    # a missing lnL must not become 0.0: lnL is negative, so 0.0 would look like a better fit
    v = float(x)
    if not math.isfinite(v):
        raise ValueError(f"non-finite numeric value: {x!r}")
    return v


def main():
    os.makedirs(OUT, exist_ok=True)
    cm = {r["HOG"]: r for r in csv.DictReader(open(CODEML_SUMMARY), delimiter="\t")}
    fl = {r["HOG"]: r for r in csv.DictReader(open(FILTER_RESULTS), delimiter="\t")}
    kept = sorted(h for h in cm if fl.get(h, {}).get("pass_all") == "1")

    rows = []
    for h in kept:
        c, q = cm[h], fl[h]
        lr_br = max(0.0, 2 * (f(c["lnL_branch"]) - f(c["lnL_M0"])))
        lr_bs = max(0.0, 2 * (f(c["lnL_bsA"]) - f(c["lnL_bsN"])))
        rows.append(dict(
            HOG=h,
            nogap_codon=c["ncodon"], aln_aa=q["aln_aa"],
            fg_dS_M0=c["fg_dS_M0"], omega_M0=c["omega_M0"],
            omega_bg=c["omega_bg"], omega_fg=c["omega_fg"],
            lnL_M0=c["lnL_M0"], lnL_branch=c["lnL_branch"],
            lnL_bsA=c["lnL_bsA"], lnL_bsN=c["lnL_bsN"],
            p2_prop=c["p2_prop"], omega2_fg=c["omega2_fg"],
            beb_n95=c["beb_n95"], beb_maxP=c["beb_maxP"],
            LRT_branch=lr_br, p_branch=chi2.sf(lr_br, 1),
            # 50:50 mixture of a point mass at 0 and chi-square with 1 df; at LRT = 0 the point
            # mass is also at or above the observed value, so the p value is 1, not 0.5
            LRT_bs=lr_bs, p_bs=(0.5 * chi2.sf(lr_bs, 1) if lr_bs > 0 else 1.0),
            accelerated=int(f(c["omega_fg"]) > f(c["omega_bg"])),
            len_ratio=q["len_ratio"],
            self_rescue=q["self_rescue"], reverse_rescue=q["reverse_rescue"]))

    # Benjamini-Hochberg, once per test over the whole set
    for key, pcol, qcol in (("branch", "p_branch", "q_branch"), ("bs", "p_bs", "q_bs")):
        qs = multipletests([r[pcol] for r in rows], alpha=FDR_Q, method="fdr_bh")[1]
        for r, qv in zip(rows, qs):
            r[qcol] = qv

    # candidate criteria, applied after the correction
    for r in rows:
        r["sig_branch"] = int(r["q_branch"] < FDR_Q and r["accelerated"])
        # branch-site also requires acceleration: model A only asks whether some codons have
        # omega2 > 1 on the foreground, so a gene slower than the background could pass
        r["sig_bs"] = int(r["q_bs"] < FDR_Q and int(r["beb_n95"] or 0) >= 1
                          and r["accelerated"])
        r["candidate"] = int(r["sig_branch"] or r["sig_bs"])

    cols = list(rows[0].keys())
    for name, sel in (("branch_results", lambda r: True),
                      ("branchsite_results", lambda r: True),
                      ("candidates", lambda r: r["candidate"])):
        with open(f"{OUT}/{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            keyf = {"branch_results": "q_branch", "branchsite_results": "q_bs",
                    "candidates": "q_bs"}[name]
            for r in sorted((x for x in rows if sel(x)), key=lambda x: x[keyf]):
                w.writerow(r)

    L = [
        "likelihood ratio tests and Benjamini-Hochberg correction",
        f"  codeml results                        {len(cm)}",
        f"  passed the filters (tested set N)     {len(kept)}",
        "",
        f"  branch      q < {FDR_Q}                   {sum(1 for r in rows if r['q_branch'] < FDR_Q)}",
        f"    of these accelerated                {sum(r['sig_branch'] for r in rows)}",
        f"  branch-site q < {FDR_Q}                   {sum(1 for r in rows if r['q_bs'] < FDR_Q)}",
        f"    of these with a BEB site P > {BEB_P}   {sum(r['sig_bs'] for r in rows)}",
        "",
        f"  candidates (branch or branch-site)    {sum(r['candidate'] for r in rows)}",
        f"    significant in both                 {sum(1 for r in rows if r['sig_branch'] and r['sig_bs'])}",
    ]
    open(f"{OUT}/summary.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
