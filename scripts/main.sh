#!/usr/bin/env bash

# # Change direction into the directory containing the scripts

# cd scripts

# # Annotate genomes

# for genome in ../assemblies/*; do
#     base=$(basename "$genome" .fasta)

#     docker run --rm -u $(id -u):$(id -g) \
#         -v $(pwd)/../assemblies:/input \
#         -v $(pwd)/../output:/output \
#         staphb/prokka \
#         prokka /input/"$(basename "$genome")" \
#             --outdir /output/ANNOTATION/prokka_"$base" \
#             --prefix "$base" \
#             --cpus 12 ## the number of cpus can be altered
# done

## Merge all cds.ffn outputs from prokka run into a single file called merged_cds.fasta
cat ../output/ANNOTATION/*/*.ffn > ../output/ANNOTATION/merged_cds.fasta

