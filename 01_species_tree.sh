#!/bin/bash
# Species tree of 12 hymenopteran genomes.
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

WORK=/PATH/TO/YOUR/WORKDIR
mkdir -p $WORK/agat $WORK/ids $WORK/primary $WORK/og $WORK/aln $WORK/trim

# species code | annotation | protein FASTA | source of the protein ID
SPECIES="
Vduc|/PATH/TO/YOUR/Vduc.braker.gff3|/PATH/TO/YOUR/Vduc.braker.aa|mrna_id
Vman|/PATH/TO/YOUR/Vman.genomic.gff|/PATH/TO/YOUR/Vman.protein.faa|cds_protein_id
Vcra|/PATH/TO/YOUR/Vcra.genomic.gff|/PATH/TO/YOUR/Vcra.protein.faa|cds_protein_id
Vvel|/PATH/TO/YOUR/Vvel.genomic.gff|/PATH/TO/YOUR/Vvel.protein.faa|cds_protein_id
Vvul|/PATH/TO/YOUR/Vvul.genomic.gff|/PATH/TO/YOUR/Vvul.protein.faa|cds_protein_id
Vpen|/PATH/TO/YOUR/Vpen.genomic.gff|/PATH/TO/YOUR/Vpen.protein.faa|cds_protein_id
Pfus|/PATH/TO/YOUR/Pfus.genomic.gff|/PATH/TO/YOUR/Pfus.protein.faa|cds_protein_id
Pdom|/PATH/TO/YOUR/Pdom.genomic.gff|/PATH/TO/YOUR/Pdom.protein.faa|cds_protein_id
Pcan|/PATH/TO/YOUR/Pcan.genomic.gff|/PATH/TO/YOUR/Pcan.protein.faa|cds_protein_id
Rbru|/PATH/TO/YOUR/Rbru.annotation.gtf|/PATH/TO/YOUR/Rbru.protein.fasta|mrna_id
Anig|/PATH/TO/YOUR/Anig.genes.gff3|/PATH/TO/YOUR/Anig.pep.fa.gz|cds_protein_id
Amel|/PATH/TO/YOUR/Amel.genomic.gff|/PATH/TO/YOUR/Amel.protein.faa|cds_protein_id
"
S12="Amel Anig Pcan Pdom Pfus Rbru Vcra Vduc Vman Vpen Vvel Vvul"

# AGAT
while IFS='|' read -r sp ann prot mode; do
  [ -z "$sp" ] && continue
  agat_sp_keep_longest_isoform.pl -gff "$ann" -o $WORK/agat/$sp.longest.gff &
done <<< "$SPECIES"
wait

# protein ID of the retained transcript
while IFS='|' read -r sp ann prot mode; do
  [ -z "$sp" ] && continue
  if [ $mode = cds_protein_id ]; then
    awk -F'\t' '$3=="CDS"{if(match($9,/protein_id=[^;]+/)) print substr($9,RSTART+11,RLENGTH-11)}' \
        $WORK/agat/$sp.longest.gff | sort -u > $WORK/ids/$sp.ids
  else
    awk -F'\t' '$3=="CDS"{if(match($9,/Parent=[^;]+/)) print substr($9,RSTART+7,RLENGTH-7)}' \
        $WORK/agat/$sp.longest.gff | sort -u > $WORK/ids/$sp.ids
  fi
done <<< "$SPECIES"

# protein extraction, species code prefixed to the header
# $WORK/primary is the OrthoFinder input directory: every file in it is read as one species,
# so the pipeline writes only the 12 final FASTAs there and keeps no intermediate files
while IFS='|' read -r sp ann prot mode; do
  [ -z "$sp" ] && continue
  if [ $sp = Anig ]; then          # this proteome carries version suffixes in the headers
    seqkit grep --id-regexp '^(\S+?)\.\d+\s' -f $WORK/ids/$sp.ids "$prot"
  else
    seqkit grep -f $WORK/ids/$sp.ids "$prot"
  fi | seqkit seq -i | seqkit replace -p '^' -r "${sp}_" > $WORK/primary/$sp.faa
done <<< "$SPECIES"

# OrthoFinder. Its gene trees vary slightly between runs on identical input, so the node labels
# and a few single-copy orthogroups shift; for the set used here, see the deposited table.
orthofinder -f $WORK/primary -t 30 -a 8 -S diamond -n Vduca12sp -o $WORK/orthofinder

# single-copy orthogroups, exactly one gene in each of the 12 species
N0=$WORK/orthofinder/Results_Vduca12sp/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
awk -F'\t' -v S="$S12" 'NR==1{n=split(S,s," "); for(i=1;i<=NF;i++) col[$i]=i; next}
     { ok=1; for(k=1;k<=n;k++){ v=$(col[s[k]]); if(v=="" || index(v,",")) {ok=0; break} }
       if(ok){ printf "%s", $1; for(k=1;k<=n;k++) printf "\t%s", $(col[s[k]]); printf "\n" } }' \
    $N0 > $WORK/single_copy_HOGs.tsv

# terminal stop codon removal, non-standard residues to X, header set to the species code
cat $WORK/primary/*.faa > $WORK/primary_all.faa
awk -F'\t' -v dir=$WORK/og '{h=$1; sub(/^N0\./,"",h); f=dir"/"h".ids"
                             for(i=2;i<=NF;i++) print $i > f; close(f)}' $WORK/single_copy_HOGs.tsv
ls $WORK/og/*.ids | xargs -P 30 -I{} bash -c '
  b=$(basename {} .ids)
  seqkit grep -f {} '"$WORK"'/primary_all.faa \
   | seqkit replace -s -p "\*$" -r "" \
   | seqkit replace -s -p "[^ACDEFGHIKLMNPQRSTVWYX]" -r X \
   | seqkit replace -p "^([A-Za-z]+)_.*" -r "\$1" > '"$WORK"'/og/$b.faa'

# alignment and trimming
ls $WORK/og/*.faa | xargs -P 30 -I{} bash -c '
  b=$(basename {} .faa)
  mafft --auto --quiet --thread 1 {} > '"$WORK"'/aln/$b.aln
  trimal -in '"$WORK"'/aln/$b.aln -out '"$WORK"'/trim/$b.fa -automated1'

# IQ-TREE. --seed defaults to the CPU clock, so the values the published run drew are given
# here; bootstrap and sCF replicates also depend on -T, which is why the thread count is fixed.
cd $WORK
iqtree -p trim -m MFP+MERGE --merge rclusterf --rclusterf 10 --alrt 1000 -B 1000 -bnni -T 30 -o Amel --seed 515998 --prefix ml
iqtree -S trim -T 30 --seed 593023 --prefix loci
iqtree -t ml.treefile --gcf loci.treefile -T 30 --seed 862601 --prefix concord_gcf
iqtree -te ml.treefile -p ml.best_scheme.nex --scfl 100 -T 30 -o Amel --seed 134372 --prefix concord_scf
