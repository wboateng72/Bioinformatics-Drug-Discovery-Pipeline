import csv
from Bio import SeqIO

# Input files
gene_presence_file = "/Users/seqafrica/Desktop/drug_discovery/latest/gene_presence_absence.csv"
clustered_proteins_file = "/Users/seqafrica/Desktop/drug_discovery/latest/clustered_proteins"
protein_fasta_file = "/Users/seqafrica/Desktop/drug_discovery/latest/merged_proteins.faa"

# Output
output_fasta_file = "/Users/seqafrica/Desktop/drug_discovery/latest/core_all_proteins.faa"

# Step 1: Identify core gene clusters (present in 100% of genomes)
print("Parsing core gene clusters...")
core_clusters = set()

with open(gene_presence_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    fieldnames = reader.fieldnames

    # Find isolate columns: everything after 'No. isolates'
    start_idx = fieldnames.index("No. isolates") + 1
    isolate_columns = fieldnames[start_idx:]

    for row in reader:
        presence_count = sum(bool(row[isolate]) for isolate in isolate_columns)
        if presence_count == len(isolate_columns):  # Present in all genomes
            core_clusters.add(row["Gene"])

print(f"Found {len(core_clusters)} core gene clusters.")

# Step 2: Parse clustered_proteins file to get protein IDs for each cluster
print("Mapping core clusters to protein IDs...")
cluster_to_proteins = {}

with open(clustered_proteins_file) as f:
    for line in f:
        if not line.strip():
            continue
        cluster, proteins = line.strip().split(":")
        protein_ids = proteins.strip().split()
        if cluster in core_clusters:
            cluster_to_proteins[cluster] = protein_ids

# Step 3: Load protein sequences from Prokka .faa
print("Indexing protein sequences...")
protein_db = SeqIO.to_dict(SeqIO.parse(protein_fasta_file, "fasta"))

# Step 4: Extract all core protein sequences (no filtering)
print("Extracting core protein sequences...")
core_seqs = []

for protein_ids in cluster_to_proteins.values():
    for pid in protein_ids:
        if pid in protein_db:
            core_seqs.append(protein_db[pid])

print(f"Writing {len(core_seqs)} core protein sequences to {output_fasta_file}")
SeqIO.write(core_seqs, output_fasta_file, "fasta")
