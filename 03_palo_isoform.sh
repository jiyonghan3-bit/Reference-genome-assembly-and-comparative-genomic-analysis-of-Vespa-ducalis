#!/bin/bash
# Isoform selection with PALO for the six Vespinae species.
# PALO: Villanueva-Canas, Laurie & Alba (2013) GBE 5(2):457-467, github.com/egenomics/palo
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

PALO=/PATH/TO/YOUR/palo                      # PALO source directory, runs on python 2.7
PY27=/PATH/TO/YOUR/python2.7
GFFREAD=/PATH/TO/YOUR/gffread
SEQKIT=/PATH/TO/YOUR/seqkit
ETEPY=/PATH/TO/YOUR/python3                  # with ete3
OF=/PATH/TO/YOUR/OrthoFinder/Results         # OrthoFinder output directory
P=/PATH/TO/YOUR/WORKDIR/palo
S6="Vduc Vman Vcra Vvel Vvul Vpen"
NCHUNK=30
mkdir -p $P/tables $P/run $P/proteomes_palo

# annotation and protein FASTA per species
ann () { case $1 in Vduc) echo /PATH/TO/YOUR/Vduc.braker.gff3;; *) echo /PATH/TO/YOUR/$1.genomic.gff;; esac; }
faa () { case $1 in Vduc) echo /PATH/TO/YOUR/Vduc.braker.aa;;   *) echo /PATH/TO/YOUR/$1.protein.faa;; esac; }

# Common ancestor node of the six species, read from the tree rather than hardcoded: node
# numbers and HOG ids are run-specific. For the orthogroups used here, see the deposited table.
NODE=$($ETEPY - "$OF/Species_Tree/SpeciesTree_rooted_node_labels.txt" $S6 <<'EOF'
import sys
from ete3 import Tree
t = Tree(sys.argv[1], format=1)
s6 = sys.argv[2:]
a = t.get_common_ancestor(s6)
assert sorted(a.get_leaf_names()) == sorted(s6), a.get_leaf_names()
print(a.name)
EOF
)

# HOGs of that node with exactly one gene in each of the six species
awk -F'\t' -v S="$S6" '
  NR==1 { n=split(S,s," "); for(i=1;i<=NF;i++) col[$i]=i
          printf "HOG\tOG"; for(k=1;k<=n;k++) printf "\t%s", s[k]; printf "\n"; next }
  { ok=1; for(k=1;k<=n;k++){ v=$(col[s[k]]); if(v=="" || index(v,",")){ ok=0; break } }
    if(ok){ printf "%s\t%s", $1, $2
            for(k=1;k<=n;k++){ v=$(col[s[k]]); sub("^" s[k] "_", "", v); printf "\t%s", v }
            printf "\n" } }' $OF/Phylogenetic_Hierarchical_Orthogroups/$NODE.tsv \
  > $P/hog_6sp_single_copy.tsv

# protein to gene table and protein length, terminal stop codon removed
for sp in $S6; do
  if [ $sp = Vduc ]; then
    $GFFREAD $(ann $sp) --table @id,@geneid -o $P/tables/$sp.gffread.tsv
    awk -F'\t' -v OFS='\t' '{print $1, $2}' $P/tables/$sp.gffread.tsv | sort -u > $P/tables/$sp.prot2gene.tsv
  else
    $GFFREAD $(ann $sp) --table @id,@geneid,protein_id -o $P/tables/$sp.gffread.tsv
    awk -F'\t' -v OFS='\t' '$3!="."{print $3, $2}' $P/tables/$sp.gffread.tsv | sort -u > $P/tables/$sp.prot2gene.tsv
  fi
  $SEQKIT replace -s -p '\*$' -r '' $(faa $sp) | $SEQKIT fx2tab -n -i -l | sort -u > $P/tables/$sp.length.tsv
done

# PALO input; short aliases because PALO truncates protein names to 18 characters
col=3
for sp in $S6; do
  awk -F'\t' -v OFS='\t' -v sp=$sp -v c=$col -v T=$P/tables '
    FILENAME==ARGV[1] { gene[$1]=$2; iso[$2]=iso[$2] " " $1; next }
    FILENAME==ARGV[2] { len[$1]=$2; next }
    FNR==1 { next }
    { g=gene[$c]; seen[g]=1; ng++; ga=sprintf("%s_g%06d", sp, ng)
      print ga, g, $1, $c > (T "/" sp ".alias_gene.tsv")
      print ga > (T "/" sp ".homolog_col")
      n=split(substr(iso[g],2), ps, " ")
      for(i=1;i<=n;i++){ np++; pa=sprintf("%s_p%06d", sp, np)
        print pa, ps[i], ga > (T "/" sp ".alias_protein.tsv")
        print ga, pa, len[ps[i]] > (T "/" sp ".species.txt") } }' \
    $P/tables/$sp.prot2gene.tsv $P/tables/$sp.length.tsv $P/hog_6sp_single_copy.tsv
  col=$((col + 1))
done
cat   $(for sp in $S6; do echo $P/tables/$sp.species.txt; done) > $P/run/species.txt
paste $(for sp in $S6; do echo $P/tables/$sp.homolog_col; done) > $P/run/homologs.txt

# PALO, one chunk per thread; rows are independent
awk -F'\t' -v R=$P/run -v N=$NCHUNK '
  { k=(NR-1)%N; d=sprintf("%s/chunk_%02d/files", R, k)
    if(!(k in made)){ system("mkdir -p " d " && ln -s ../../species.txt " d "/species.txt"); made[k]=1 }
    print > (d "/homologs.txt") }' $P/run/homologs.txt
ls -d $P/run/chunk_* | xargs -P $NCHUNK -I{} bash -c \
  'cd {} && OMP_NUM_THREADS=1 '"$PY27"' '"$PALO"'/palo2.py'
cat $P/run/chunk_*/files/palo_combs.txt        > $P/run/palo_combs.txt
cat $P/run/chunk_*/files/unprocessed_genes.txt > $P/run/unprocessed_genes.txt

# aliases back to the original ids
awk -F'\t' -v OFS='\t' -v S="$S6" -v U=$P/palo_unprocessed.tsv '
  BEGIN { n=split(S,s," ") }
  FILENAME ~ /alias_protein\.tsv$/ { prot[$1]=$2; pgene[$1]=$3; next }
  FILENAME ~ /alias_gene\.tsv$/    { hog[$1]=$3; next }
  FILENAME ~ /hog_6sp_single_copy\.tsv$/ { if(FNR>1) og[$1]=$2; next }
  FILENAME ~ /homologs\.txt$/      { for(k=1;k<=n;k++) row[$1,k]=$k; next }
  FILENAME ~ /unprocessed_genes\.txt$/ { m=split($0, u, " "); r=u[3]; for(i=4;i<=m;i++) r=r " " u[i]
                                          print hog[u[1]], u[2], r > U; next }
  { if(FNR==1 && !hdr){ printf "HOG\tOG"; for(k=1;k<=n;k++) printf "\t%s", s[k]; printf "\n"; hdr=1 }
    h=hog[$1]; printf "%s\t%s", h, og[h]
    for(k=1;k<=n;k++) printf "\t%s", prot[$(k+1)]
    printf "\n" }' \
  $P/tables/*.alias_protein.tsv $P/tables/*.alias_gene.tsv $P/hog_6sp_single_copy.tsv \
  $P/run/homologs.txt $P/run/unprocessed_genes.txt $P/run/palo_combs.txt \
  > $P/palo_selected_isoforms.tsv

# selected proteins
col=3
for sp in $S6; do
  awk -F'\t' -v c=$col 'NR>1{print $c}' $P/palo_selected_isoforms.tsv > $P/proteomes_palo/$sp.ids
  $SEQKIT grep -f $P/proteomes_palo/$sp.ids $(faa $sp) | $SEQKIT seq -i \
    | $SEQKIT replace -p '^' -r "${sp}_" > $P/proteomes_palo/$sp.faa
  col=$((col + 1))
done
