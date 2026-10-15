#!/bin/bash
# CDS of the selected isoforms, kept only when the translation matches the protein.
# Replace every /PATH/TO/YOUR/... with your own path before running.
set -euo pipefail

PY=/PATH/TO/YOUR/python3                     # with Biopython
CDS_DIR=/PATH/TO/YOUR/cds                    # <species>.cds.fna, RefSeq cds_from_genomic or braker.codingseq
PROT_DIR=/PATH/TO/YOUR/proteomes             # <species>.faa
SELECTED=/PATH/TO/YOUR/WORKDIR/palo/palo_selected_isoforms.tsv
O=/PATH/TO/YOUR/WORKDIR/cds
mkdir -p $O/hog

$PY - $CDS_DIR $PROT_DIR $SELECTED $O <<'EOF'
import sys, re
from Bio.Seq import Seq
CDS_DIR, PROT_DIR, SELECTED, O = sys.argv[1:]
S6 = ["Vduc", "Vman", "Vcra", "Vvel", "Vvul", "Vpen"]
STOP = ("TAA", "TAG", "TGA")

def read_fasta(path, key=None):
    d, h, b = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None: d[h] = "".join(b)
            h = key(line) if key else line[1:].split()[0]; b = []
        else:
            b.append(line.strip())
    if h is not None: d[h] = "".join(b)
    return d

# RefSeq CDS headers carry [protein_id=...]; BRAKER coding sequences use the transcript id
pid_re = re.compile(r"\[protein_id=([^\]]+)\]")
cds, prot = {}, {}
for sp in S6:
    cds[sp] = read_fasta(f"{CDS_DIR}/{sp}.cds.fna",
                         key=lambda l: (pid_re.search(l).group(1) if "[protein_id=" in l
                                        else l[1:].split()[0]))
    prot[sp] = read_fasta(f"{PROT_DIR}/{sp}.faa")

rows = [l.rstrip("\n").split("\t") for l in open(SELECTED)]
hdr, rows = rows[0], rows[1:]
col = {sp: hdr.index(sp) for sp in S6}

excl = []
for r in rows:
    hog = r[0]; seqs = {}; bad = None
    for sp in S6:
        pid = r[col[sp]]
        p = prot[sp][pid].rstrip("*")
        s = cds[sp].get(pid, "").upper()
        if not s:          bad = (sp, pid, "no CDS"); break
        if s[-3:] in STOP: s = s[:-3]
        if len(s) % 3:     bad = (sp, pid, "length not a multiple of three"); break
        t = str(Seq(s).translate(table=1))
        if "*" in t:       bad = (sp, pid, "internal stop codon"); break
        if t != p:         bad = (sp, pid, "translation does not match the protein"); break
        seqs[sp] = s
    if bad:
        excl.append((hog, *bad)); continue
    with open(f"{O}/hog/{hog}.cds.fa", "w") as f:
        for sp in S6:
            f.write(">%s\n%s\n" % (sp, seqs[sp]))

with open(f"{O}/excluded.tsv", "w") as f:
    f.write("HOG\tspecies\tprotein_id\treason\n")
    for e in excl:
        f.write("\t".join(e) + "\n")
EOF
