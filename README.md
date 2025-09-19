
1. **Annotation of genomes**



The assembled genomes were annotated using prokka using the command : 
prokka klebsiella.fasta \
  --outdir prokka_klebsiella \
  --prefix klebsiella \
  --genus Klebsiella --species pneumoniae \
  --locustag KPN --usegenus --cpus 12







The protein.faa outputs were assembled in a single file called merged_proteins using the command: cat */*.faa > merged_proteins.faa





The cds.ffn outputs were assembled in a single file called merged_cds using the command: cat */*.faa > merged_cds.fa




The pan genome was generated using Roary using the command:
roary -e -n -v -p 12 -i 90 -cd 100 -f pangenome *.gff



The CD-HIT program was used to cluster the proteins based on similarity using a threshold of 0.7 and word size 5 using the command: cd-hit -i core_proteins.faa -o clustered_core_proteins.faa -c 0.7 -n 5



The CD-HIT program was used to cluster the cds based on similarity using a threshold of 0.8 and word size 5 using the command: cd-hit-est -i core_cds.fa -o clustered_core_cds.fa -c 0.8 -n 5



The VFDB database in Abricate was used to predict the virulence genes composition of the genomes using the clustered_core_cds.fa. The command used was: abricate --minid 80 --threads 12 clustered_core_cds.fa > virulence_genes.csv

The human genome was downloaded using the command: wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz

It was gunzipped



The database was created using the command: makeblastdb -in UP000005640_9606.fasta -dbtype prot -out human9606_db

Blasting was done using the command: blastp -query clustered_core_proteins.faa -db human9606_db -out results.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5 -num_threads 12

Sequences homologous to the human genome was filtered out (threshold- >= 30% identity, >=100 bitscore and <1e-3)




The bacteria protein sequences in the essential database was downloaded (http://origin.tubic.org/deg/public/index.php)

The database was created using the command: makeblastdb -in DEG10.aa.fa -dbtype prot -out DEG_db

The command below was used to blast the filtered_human_homolog protein file:

(echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"; blastp -query human_homolog_filtered_proteins.fa -db DEG_db -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5 -num_threads 12) DEG_results.txt

Interproscan was used to assess the protein families,  domains and motifs of the essential genes.

KAAS was used to assign KO numbers to the proteins.

Psortb was used to determine the localisation of the sequences (https://psort.org/psortb/)
