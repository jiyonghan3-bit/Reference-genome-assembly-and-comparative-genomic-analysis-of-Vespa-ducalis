#!/bin/bash
# Codon alignment with PRANK and gap column removal with trimAl.
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

PRANK=/PATH/TO/YOUR/prank
TRIMAL=/PATH/TO/YOUR/trimal
ETEPY=/PATH/TO/YOUR/python3                  # with ete3
ML_TREEFILE=/PATH/TO/YOUR/ml.treefile        # 12 species ML tree
CDS=/PATH/TO/YOUR/WORKDIR/cds/hog            # <orthogroup>.cds.fa
O=/PATH/TO/YOUR/WORKDIR/alignment
mkdir -p $O/prank $O/nogap $O/tmp

# guide tree, the six species pruned from the species tree with branch lengths kept
$ETEPY - $ML_TREEFILE $O/guide_6sp.nwk <<'EOF'
import sys
from ete3 import Tree
ml, out = sys.argv[1:]
S6 = ["Vduc", "Vman", "Vcra", "Vvel", "Vvul", "Vpen"]
t = Tree(ml, format=1)
t.prune(S6, preserve_branch_length=True)
assert sorted(t.get_leaf_names()) == sorted(S6), t.get_leaf_names()
t.write(outfile=out, format=5)
EOF

# alignment and trimming; PRANK writes into the working directory, so each gene runs in its own
ls $CDS/*.cds.fa | xargs -P 40 -I{} bash -c '
  f={}; b=$(basename $f .cds.fa); w='"$O"'/tmp/$b
  mkdir -p $w
  ( cd $w && '"$PRANK"' -d=$f -t='"$O"'/guide_6sp.nwk -o=$b -codon -F -once )
  mv $w/$b.best.fas '"$O"'/prank/$b.best.fas
  '"$TRIMAL"' -in '"$O"'/prank/$b.best.fas -out '"$O"'/nogap/$b.fa -nogaps'
# the per-gene directories under $O/tmp hold PRANK's other output files and are left in place
