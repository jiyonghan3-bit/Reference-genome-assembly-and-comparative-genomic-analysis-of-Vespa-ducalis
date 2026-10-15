#!/bin/bash
# codeml, four models per orthogroup: M0, two-ratio branch, branch-site A and its null.
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

PY=/PATH/TO/YOUR/python3                     # with Biopython
ETEPY=/PATH/TO/YOUR/python3                  # with ete3
CODEML=/PATH/TO/YOUR/codeml
WORKER=/PATH/TO/YOUR/06_codeml_worker.py
ML_TREEFILE=/PATH/TO/YOUR/ml.treefile        # 12 species ML tree
A=/PATH/TO/YOUR/WORKDIR/alignment/nogap      # codon alignments, gap columns removed
O=/PATH/TO/YOUR/WORKDIR/codeml
mkdir -p $O

# unrooted six species tree with the foreground branch marked; the " 6 1" header line is required
$ETEPY - $ML_TREEFILE $O/tree_6sp_fg.nwk <<'EOF'
import sys
from ete3 import Tree
ml, out = sys.argv[1:]
S6 = ["Vduc", "Vman", "Vcra", "Vvel", "Vvul", "Vpen"]
t = Tree(ml, format=1)
t.prune(S6)
t.unroot()
assert len(t.children) == 3, t.children
assert sorted(t.get_leaf_names()) == sorted(S6), t.get_leaf_names()
nwk = t.write(format=9).replace("Vduc", "Vduc #1")
assert nwk.count("#1") == 1, nwk
open(out, "w").write(" 6 1\n" + nwk + "\n")
EOF

# one worker per orthogroup
export ALNDIR=$A TREE=$O/tree_6sp_fg.nwk OUTDIR=$O CODEML=$CODEML
ls $A/*.fa | xargs -n1 basename | sed 's/\.fa$//' | sort \
  | xargs -P 50 -I{} $PY $WORKER {} \
    > $O/codeml_results.tsv 2> $O/failed.tsv

# column names; 08_lrt_fdr.py reads codeml_results.header.tsv
HDR="HOG\tncodon\tfg_branch\tlnL_M0\tomega_M0\tfg_dS_M0\tfg_dN_M0\tfg_t_M0\tlnL_branch\tomega_bg\tomega_fg\tlnL_bsA\tlnL_bsN\tp2_prop\tomega2_fg\tbeb_n\tbeb_n95\tbeb_maxP\ttrials_M0\ttrials_branch\ttrials_bsA\ttrials_bsN\tsec\tstatus"
{ printf "$HDR\n"; cat $O/codeml_results.tsv; } > $O/codeml_results.header.tsv
cat $O/*/beb.tsv 2>/dev/null | sort -k1,1 -k2,2n > $O/beb_all.tsv
