# The pipeline generally consists of:

- **Annotation of genomes**

- **Pangenome analysis**

- **Selection of core genes and proteins**

- **Filtering of human homologs**

- **Pathway/metabolic analysis**

- **Subcellular localization prediction**

- **Virulence/AMR analysis**

- **Druggability scoring**


# Installation of tools

### The NCBI AMRfinder plus docker was set up using the command
```
docker pull ncbi/amr
```

# Installation of databases

### The human genome was downloaded and set up using the commands
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz

makeblastdb -in UP000005640_9606.fasta -dbtype prot -out human9606_db
```













Blasting was done using the command: blastp -query clustered_core_proteins.faa -db human9606_db -out results.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5 -num_threads 12

Sequences homologous to the human genome was filtered out (threshold- >= 30% identity, >=100 bitscore and <1e-3)




The bacteria protein sequences in the essential database was downloaded (http://origin.tubic.org/deg/public/index.php)

The database was created using the command: makeblastdb -in DEG10.aa.fa -dbtype prot -out DEG_db

The command below was used to blast the filtered_human_homolog protein file:

(echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"; blastp -query human_homolog_filtered_proteins.fa -db DEG_db -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5 -num_threads 12) DEG_results.txt

Interproscan was used to assess the protein families,  domains and motifs of the essential genes.

KAAS was used to assign KO numbers to the proteins.

Psortb was used to determine the localisation of the sequences (https://psort.org/psortb/)
