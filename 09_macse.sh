#!/bin/bash
# Alignment sensitivity: the candidates are realigned with MACSE and tested again.
# Significance uses the raw p thresholds of the main analysis, not a new BH correction.
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

MACSE=/PATH/TO/YOUR/macse
TRIMAL=/PATH/TO/YOUR/trimal
PY=/PATH/TO/YOUR/python3                     # with Biopython and scipy
CODEML=/PATH/TO/YOUR/codeml
WORKER=/PATH/TO/YOUR/06_codeml_worker.py
CDS=/PATH/TO/YOUR/WORKDIR/cds/hog            # <orthogroup>.cds.fa
TREE=/PATH/TO/YOUR/WORKDIR/codeml/tree_6sp_fg.nwk
TESTS=/PATH/TO/YOUR/WORKDIR/tests            # branch_results.tsv, candidates.tsv
O=/PATH/TO/YOUR/WORKDIR/macse
mkdir -p $O/aln $O/codeml $O/tmp

# candidates
awk -F'\t' 'NR>1{print $1}' $TESTS/candidates.tsv | sort > $O/candidates.txt

# MACSE realignment; frameshifts and internal stops are replaced so that codeml can read them
xargs -a $O/candidates.txt -P 10 -I{} bash -c '
  h={}; T='"$O"'/tmp/$h
  mkdir -p $T
  '"$MACSE"' -prog alignSequences -seq '"$CDS"'/$h.cds.fa -gc_def 1 \
      -out_NT $T/raw_NT.fa -out_AA $T/raw_AA.fa > $T/macse.log 2>&1
  '"$MACSE"' -prog exportAlignment -align $T/raw_NT.fa \
      -codonForInternalStop NNN -codonForFinalStop --- -codonForInternalFS --- -codonForExternalFS --- \
      -out_NT $T/clean_NT.fa -out_AA $T/clean_AA.fa >> $T/macse.log 2>&1
  '"$TRIMAL"' -in $T/clean_NT.fa -out '"$O"'/aln/$h.fa -nogaps 2>/dev/null'

# codeml, same worker and same tree as the main analysis
export ALNDIR=$O/aln TREE=$TREE OUTDIR=$O/codeml CODEML=$CODEML
ls $O/aln/*.fa | xargs -n1 basename | sed 's/\.fa$//' \
  | xargs -P 10 -I{} $PY $WORKER {} > $O/codeml_results.tsv 2> $O/failed.tsv

# comparison with the main analysis
$PY - $TESTS/branch_results.tsv $TESTS/candidates.tsv $O/codeml_results.tsv $O <<'EOF'
import sys, csv
from scipy.stats import chi2
allres, cand, macse, out = sys.argv[1:]

# largest p that was still q < 0.05 in the main analysis
pb = ps = 0.0
for r in csv.DictReader(open(allres), delimiter="\t"):
    if float(r["q_branch"]) < 0.05: pb = max(pb, float(r["p_branch"]))
    if float(r["q_bs"]) < 0.05:     ps = max(ps, float(r["p_bs"]))

C = {r["HOG"]: r for r in csv.DictReader(open(cand), delimiter="\t")}
COL = "HOG ncodon fg_branch lnL_M0 omega_M0 fg_dS_M0 fg_dN_M0 fg_t_M0 lnL_branch omega_bg omega_fg \
lnL_bsA lnL_bsN p2_prop omega2_fg beb_n beb_n95 beb_maxP \
trials_M0 trials_branch trials_bsA trials_bsN sec status".split()
M = {}
for line in open(macse):
    v = line.rstrip("\n").split("\t")
    if len(v) == len(COL) and v[-1] == "OK": M[v[0]] = dict(zip(COL, v))

rows = []
for h, c in C.items():
    m = M.get(h)
    if not m:
        rows.append(dict(HOG=h, macse="failed")); continue
    lb = max(0.0, 2*(float(m["lnL_branch"]) - float(m["lnL_M0"])))
    ls = max(0.0, 2*(float(m["lnL_bsA"]) - float(m["lnL_bsN"])))
    # same mixture handling as 08: at LRT = 0 the branch-site p value is 1, not 0.5
    p_b = chi2.sf(lb, 1)
    p_s = 0.5 * chi2.sf(ls, 1) if ls > 0 else 1.0
    acc = float(m["omega_fg"]) > float(m["omega_bg"])
    sb = int(p_b <= pb and acc)
    ss = int(p_s <= ps and int(m["beb_n95"]) >= 1)
    rows.append(dict(HOG=h, macse="OK",
        prank_codon=c["nogap_codon"], macse_codon=m["ncodon"],
        prank_omega_fg=c["omega_fg"], macse_omega_fg=m["omega_fg"],
        prank_q_branch=c["q_branch"], macse_p_branch=f"{p_b:.4g}",
        prank_q_bs=c["q_bs"], macse_p_bs=f"{p_s:.4g}",
        prank_beb=c["beb_n95"], macse_beb=m["beb_n95"],
        prank_sig=("branch" if c["sig_branch"]=="1" else "")+("bs" if c["sig_bs"]=="1" else ""),
        macse_sig=("branch" if sb else "")+("bs" if ss else ""),
        # robust means the SAME test is significant under both alignments
        robust=int((c["sig_branch"]=="1" and sb) or (c["sig_bs"]=="1" and ss))))
cols = list(max(rows, key=len).keys())
with open(f"{out}/macse_comparison.tsv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader()
    for r in sorted(rows, key=lambda x: -x.get("robust",0)): w.writerow(r)
n_ok = sum(1 for r in rows if r.get("macse")=="OK")
n_rb = sum(r.get("robust",0) for r in rows)
L = ["MACSE alignment sensitivity",
     f"  candidates {len(rows)}, realigned and tested {n_ok}",
     f"  supported by the same test under both alignments {n_rb}",
     f"  not supported by the same test under both alignments {n_ok - n_rb}"]
open(f"{out}/summary.txt","w").write("\n".join(L)+"\n")
print("\n".join(L))
EOF
# the per-gene directories under $O/tmp hold the MACSE intermediates and are left in place
