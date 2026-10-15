# Vespa ducalis: phylogenomics, divergence times and lineage-specific rate acceleration

Scripts for the phylogenomic and dN/dS analyses reported in the accompanying manuscript. These document the study workflow rather than a fully automated pipeline. Replace `/PATH/TO/YOUR/...` with local paths before running.

## Workflow

Run 01 first, followed by 02 and 03–09. Both 02 and the 03–09 chain use outputs from 01; 03–09 does not depend on 02. Helper scripts are called by the numbered scripts.

1. Select the longest protein isoform per gene for 12 species, infer orthogroups, align single-copy genes, and infer the ML species tree with gene and site concordance factors.
2. Select clock-like genes, concatenate alignments, calculate an LG+G4 Hessian, and estimate divergence times with fossil calibrations and sensitivity analyses.
3. Extract single-copy HOGs for the six Vespinae species and select isoform combinations that minimise length differences.
4. Extract CDS and retain sequences whose translations match the selected proteins.
5. Generate codon alignments using a pruned guide tree and filter gap-containing positions while preserving the reading frame.
6. Specify the foreground branch and fit four codeml models per gene group using multiple starting values.
7. Apply alignment and isoform integrity filters.
8. Perform likelihood ratio tests, apply Benjamini–Hochberg correction, and identify candidates.
9. Realign candidates with MACSE and repeat the codeml analyses to assess sensitivity to alignment method.

Step 07 uses alignments and annotations, without using test statistics or significance values. Retained groups define the multiple-testing set in 08. Candidate criteria are applied after correction; 09 is a follow-up analysis of candidates rather than a second genome-wide screen.

The public scripts omit the original resume logic and retain intermediate files. Use fresh output directories for new runs. The OrthoFinder input directory in 01 must contain only the 12 intended protein FASTAs. Figure code, assembly, repeat and structural annotation, and symbiont-genome analyses are outside this repository's scope.

## Reproducibility

Three OrthoFinder 2.5.5 runs with identical inputs produced 6,204, 6,205 and 6,206 single-copy HOGs at the root. Comparisons identified at least 6,201 groups with identical protein membership, and the species-tree topology was unchanged. DIAMOND results and MCL clustering were identical. Downstream differences were traced to Python set iteration order, which changed sequence order in DendroBLAST distance matrices and affected gene-tree inference where distances were tied.

Two further runs with `PYTHONHASHSEED=0` produced identical distance matrices and gene trees in the tested environment. Set this variable before launching OrthoFinder for new runs. This controls the variation identified here, but does not guarantee identical results across software versions or environments.

The manuscript analysis did not fix the hash seed. Use the accompanying orthogroup tables to repeat downstream analyses with the original gene assignments; fixing the seed now does not recover the original run. OG identifiers were unchanged in our comparisons, whereas HOG identifiers and internal node labels are run-specific. Compare groups by protein membership and nodes by descendant taxa.

The original IQ-TREE seeds and thread settings are recorded in 01. Retain these settings when repeating the analyses.

## Software and parallel execution

AGAT 0.8.1, SeqKit 2.3.0, OrthoFinder 2.5.5, DIAMOND 2.1.24, MAFFT 7.525, trimAl 1.5.1, IQ-TREE 3.0.1, SortaDate commit `15bbbc4`, FASconCAT-G 1.06.1, PALO commit `008638fb80ca`, gffread 0.12.7, PRANK v.170427, PAML 4.10.10, MACSE 2.07, Biopython 1.86, ETE3 3.1.3, SciPy 1.15.2 and statsmodels 0.14.6.

The AGAT 0.8.1 installation reports `v0.8.0` in its banner. Unmodified PALO runs in a separate environment with Python 2.7.15 and NumPy 1.16.5; the remaining Python scripts use Python 3. External dependencies include MCL and FastME for OrthoFinder and Java for MACSE.

Original parallelisation settings: OrthoFinder `-t 30 -a 8`; MAFFT/trimAl, 30 concurrent jobs with MAFFT `--thread 1`; PALO, 30 chunks and 30 concurrent jobs with `OMP_NUM_THREADS=1`; PRANK, 40 jobs; codeml in 06, 50 jobs; step 07, worker count supplied as an argument (default 40); MACSE/codeml in 09, 10 jobs. Adjust to available resources, or retain the original settings for reproducibility comparisons.

## Input data

Comparison genomes and annotations were obtained from RefSeq, Ensembl and Figshare under the accessions listed in the manuscript. Decompress the Ensembl GFF3 for A. nigricornis before running 01. Accessions for the V. ducalis assembly, annotation and reads are also provided in the manuscript.

The accompanying orthogroup tables record the assignments used for 01 and 03–09. The root-level HOG table and the table for the six-species Vespinae clade represent different extraction levels and should not be substituted for one another.
