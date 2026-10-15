#!/bin/bash
# Divergence time estimation for the 12 hymenopteran species.
# Requires the gene trees, the species tree and the trimmed alignments from 01_species_tree.sh.
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

WORK=/PATH/TO/YOUR/WORKDIR
TRIM=/PATH/TO/YOUR/trim                                  # trimmed single-gene alignments
LOCI_TREEFILE=/PATH/TO/YOUR/loci.treefile                # gene trees
LOCI_SCHEME=/PATH/TO/YOUR/loci.best_scheme.nex           # gene order of loci.treefile
ML_TREEFILE=/PATH/TO/YOUR/ml.treefile                    # species tree
SORTADATE=/PATH/TO/YOUR/SortaDate/src
FASCONCAT=/PATH/TO/YOUR/FASconCAT-G_v1.06.1.pl
LG_DAT=/PATH/TO/YOUR/paml/dat/lg.dat
mkdir -p $WORK/01_gene_selection/genetrees $WORK/02_matrix/genes $WORK/03_mcmctree/hessian

# phyx: one file per gene tree, rooted on Apis mellifera
O=$WORK/01_gene_selection
grep -o 'charset [^ ]*' $LOCI_SCHEME | sed 's/charset //; s/\.fa$//' > $O/gene_order.txt
paste $O/gene_order.txt $LOCI_TREEFILE | while IFS=$'\t' read -r g tr; do
  printf "%s\n" "$tr" > $O/unrooted.tmp
  pxrr -t $O/unrooted.tmp -g Amel > $O/genetrees/$g.tre
done
rm -f $O/unrooted.tmp

# SortaDate: root-to-tip variance, bipartition support against the species tree, top 500 genes
python $SORTADATE/get_var_length.py $O/genetrees --flend .tre --outg Amel --outf $O/var.txt
python $SORTADATE/get_bp_genetrees.py $O/genetrees $ML_TREEFILE --flend .tre --outf $O/bp.txt
python $SORTADATE/combine_results.py $O/var.txt $O/bp.txt --outf $O/comb.txt
python $SORTADATE/get_good_genes.py $O/comb.txt --max 500 --outf $O/good_genes.txt
awk 'NR>1{sub(/\.tre$/,"",$1); print $1}' $O/good_genes.txt > $O/selected_500.txt

# FASconCAT-G: concatenation of the selected alignments
M=$WORK/02_matrix
while read -r g; do cp $TRIM/$g.fa $M/genes/$g.fas; done < $O/selected_500.txt
( cd $M/genes && perl $FASCONCAT -s -p -p -l )
mv $M/genes/FcC_* $M/
awk 'NR==1{print " "$1" "$2; next} {printf "%s  %s\n", $1, $2}' $M/FcC_supermatrix.phy > $M/mcmctree_aa.phy

# calibrated trees: main analysis, Palaeovespa alternative, Symmorphus crown alternative
C=$WORK/03_mcmctree
SEQ=$M/mcmctree_aa.phy
cat > $C/main.trees <<'EOF'
12 1
(Amel, ((Anig, Rbru), (((Pcan, Pfus), Pdom), (((Vcra, Vvel), (Vduc, Vman)), (Vpen, Vvul))) '>0.60') '>0.90') 'U(1.576, 0.05)';
EOF
cat > $C/alt.trees <<'EOF'
12 1
(Amel, ((Anig, Rbru), (((Pcan, Pfus), Pdom), (((Vcra, Vvel), (Vduc, Vman)), (Vpen, Vvul)) '>0.60')) '>0.90') 'U(1.576, 0.05)';
EOF
cat > $C/eum.trees <<'EOF'
12 1
(Amel, ((Anig, Rbru) '>0.90', (((Pcan, Pfus), Pdom), (((Vcra, Vvel), (Vduc, Vman)), (Vpen, Vvul))) '>0.60')) 'U(1.576, 0.05)';
EOF

ctl_common () {   # $1=treefile  $2=usedata  $3=clock  $4=seed
cat <<EOF
          seed = $4
       seqfile = $SEQ
      treefile = $1
      mcmcfile = mcmc.txt
       outfile = out.txt
         ndata = 1
       seqtype = 2
       usedata = $2
         clock = $3
       BDparas = 1 1 0.1 M
   rgene_gamma = 2 9 1
  sigma2_gamma = 1 10 1
      finetune = 1: .1 .1 .1 .1 .1 .1
         print = 1
        burnin = 100000
      sampfreq = 200
       nsample = 20000
EOF
}

# MCMCTree usedata = 3, then codeml with LG+G4 for the branch lengths, gradient and Hessian
cd $C/hessian
ctl_common $C/main.trees 3 2 -1 > mcmctree.ctl
# usedata = 3 makes mcmctree write tmp0001.* and then invoke codeml itself. That call is not relied
# on, so its exit status is ignored; the Hessian that in.BV uses is recomputed with LG+G4 below.
mcmctree mcmctree.ctl || true
[ -s tmp0001.ctl ] || { echo "mcmctree did not write tmp0001.ctl"; exit 1; }

# set the substitution model: delete these lines if present, then append, so they are set either way
sed -i -E '/^ *(model|aaRatefile|fix_alpha|alpha|ncatG) *=/d' tmp0001.ctl
printf ' model = 2\n aaRatefile = %s\n fix_alpha = 0\n alpha = 0.5\n ncatG = 4\n' "$LG_DAT" >> tmp0001.ctl
codeml tmp0001.ctl

grep -q "dGamma (ncatG=4)" tmp0001.out || { echo "LG+G4 was not applied"; exit 1; }
[ -s rst2 ] || { echo "codeml produced no rst2"; exit 1; }
cp rst2 $C/in.BV      # rst2, not out.BV, which may hold a default-model Hessian

# MCMC: independent rates x2, autocorrelated rates x2, prior only, two alternative calibrations
run () {   # $1=name  $2=treefile  $3=usedata  $4=clock  $5=seed
  local d=$C/$1; mkdir -p $d
  ( cd $d && ln -sf $C/in.BV in.BV && ctl_common $2 "$3" $4 $5 > mcmctree.ctl && mcmctree mcmctree.ctl ) &
  PIDS+=($!)
}
PIDS=()
run ir_run1 $C/main.trees "2 in.BV" 2 -1
run ir_run2 $C/main.trees "2 in.BV" 2 12345
run ar_run1 $C/main.trees "2 in.BV" 3 -1
run ar_run2 $C/main.trees "2 in.BV" 3 54321
run prior   $C/main.trees 0          2 -1
run alt_ir  $C/alt.trees  "2 in.BV"  2 -1
run alt_eum $C/eum.trees  "2 in.BV"  2 -1
wait "${PIDS[@]}"

# MCMC: birth-death time prior with the conditional (C) scheme
mkdir -p $C/ir_bdC
cd $C/ir_bdC
ln -sf $C/in.BV in.BV
sed 's|BDparas = 1 1 0.1 M|BDparas = 1 1 0.1 C|' $C/ir_run1/mcmctree.ctl > mcmctree.ctl
mcmctree mcmctree.ctl
