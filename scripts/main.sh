#!/usr/bin/env bash

#cd ..

# # Annotate genomes

# for genome in assemblies/*; do
#     base=$(basename "$genome" .fasta)

#     docker run --rm -u $(id -u):$(id -g) \
#         -v "$(pwd)/assemblies":/input \
#         -v "$(pwd)/output":/output \
#         staphb/prokka \
#         prokka /input/"$(basename "$genome")" \
#             --outdir /output/ANNOTATION/prokka_"$base" \
#             --prefix "$base" \
#             --cpus 12 ## the number of cpus can be altered
# done

### Pangenome analysis
# mkdir -p output/ROARY/GFF
# cp output/ANNOTATION/*/*.gff output/ROARY/GFF

# #Run Roary pangenome analysis
# -e -n   : fast core gene alignment with MAFFT
# -v      : verbose mode
# -p 12   : use 12 threads
# -i 95   : minimum blastp identity 95% for orthologous clusters
# -cd 95  : genes or protein found in 95% of isolates
# -f /output/RESULTS : output folder
# /input/*.gff : input files

# docker run --rm -u "$(id -u):$(id -g)" \
#   -v "$(realpath output/ROARY/GFF)":/input \
#   -v "$(realpath output/ROARY)":/output \
#   sangerpathogens/roary \
#   bash -lc "roary -e -n -v -p 12 -i 95 -cd 95 -f /output/RESULTS /input/*.gff"

###Create core proteins and cds files 

## Merge all cds.ffn outputs from prokka run into a single file called merged_cds.faa
# mkdir output/CORE_PROT_CDS
# cat output/ANNOTATION/*/*.faa > output/CORE_PROT_CDS/merged_proteins.faa
# cat output/ANNOTATION/*/*.ffn > output/CORE_PROT_CDS/merged_cds.fasta

## Cluster the protein files (.faa) files using the clustering threshold set when using Roary
# docker run --rm -u "$(id -u):$(id -g)" \
#   -v "$(realpath output/ROARY/GFF)":/input \
#   -v "$(realpath output/ROARY)":/output \
#   sangerpathogens/roary \
#   bash -lc "cd-hit -i output/CORE_PROT_CDS/merged_proteins.faa -o output/CORE_PROT_CDS/clustered_proteins.faa -c 0.95 -n 4" #the right n (word size) helps prevent redunduncy

## Copy all files needed for core protein selection into CORE_PROT_CDS folder
#cp output/ROARY/RESULTS/clustered_proteins output/ROARY/RESULTS/gene_presence_absence.csv output/CORE_PROT_CDS/

##Run core_proteins.py script to generate file (core_proteins.faa) containing core protein sequences (edit the script to define core genes)
#python3 scripts/core_proteins.py

## Copy the headers of the core genes in core_proteins.faa into a file called core_protein_headers.csv 
#grep "^>" output/CORE_PROT_CDS/core_proteins.faa > output/CORE_PROT_CDS/core_protein_headers.csv

##Run cds_protein_id_matches.py script to generate file (core_cds.fasta) containing core cds sequences
#python3 scripts/cds_protein_id_matches.py
