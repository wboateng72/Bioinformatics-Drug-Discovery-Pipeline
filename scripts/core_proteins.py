import csv
from Bio import SeqIO

# Input files
gene_presence_file = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/gene_presence_absence.csv"
clustered_proteins_file = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/clustered_proteins"
clustered_protein_fasta = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/clustered_proteins.faa"
output_fasta = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/core_proteins.faa"

import csv
from Bio import SeqIO

# Input files
gene_presence_file = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/gene_presence_absence.csv"
clustered_proteins_file = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/clustered_proteins"
clustered_protein_fasta = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/clustered_proteins.faa"
output_fasta = "/Users/seqafrica/PERSONAL_WORKSPACE/NOGUCHI FILES/GITHUB REPOSITORY/Bioinformatics-Drug-Discovery-Pipeline/output/CORE_PROT_CDS/core_proteins.faa"

# Step 1. Count isolates = number of columns after "Avg group size nuc"
with open(gene_presence_file, newline="") as f:
    reader = csv.DictReader(f)
    headers = reader.fieldnames
    avg_index = headers.index("Avg group size nuc")
    isolate_columns = headers[avg_index + 1:]
    total_isolates = len(isolate_columns)
    print(f"DEBUG: Found {total_isolates} isolates in dataset.")

# Step 2. Identify core genes (>=95%)
core_genes = set()
with open(gene_presence_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            no_isolates = int(row["No. isolates"])
        except ValueError:
            continue
        percentage = (no_isolates / total_isolates) * 100
        if percentage >= 95:
            core_genes.add(row["Gene"])

print(f"DEBUG: Found {len(core_genes)} core genes.")

# Step 3. Match core genes in clustered_proteins and take first representative ID
core_reps = set()
with open(clustered_proteins_file) as f:
    for line in f:
        if not line.strip():
            continue
        gene, proteins = line.split(":")
        gene = gene.strip()
        if gene in core_genes:
            rep_id = proteins.strip().split()[0]  # take the first protein after colon
            core_reps.add(rep_id)

print(f"DEBUG: Extracted {len(core_reps)} representative IDs.")

# Step 4. Extract sequences whose headers contain any of the representative IDs
records = []
for record in SeqIO.parse(clustered_protein_fasta, "fasta"):
    for rep in core_reps:
        if rep in record.description:  # check substring match
            records.append(record)
            break

SeqIO.write(records, output_fasta, "fasta")
print(f"Saved {len(records)} core protein sequences to {output_fasta}")
