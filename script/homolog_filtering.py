import sys

def filter_sequences(input_fasta_file, ids_to_exclude_file, output_file):
    """
    Filters sequences from a FASTA file, excluding those whose headers
    partially match IDs from a given CSV file.

    Args:
        input_fasta_file (str): Path to the input FASTA file.
        ids_to_exclude_file (str): Path to the CSV file containing IDs to exclude.
        output_file (str): Path to the output file for filtered sequences.
    """
    try:
        # Read IDs to exclude from the CSV file
        with open(ids_to_exclude_file, 'r') as f:
            ids_to_exclude = {line.strip() for line in f if line.strip()}

        # Read and filter sequences from the FASTA file
        with open(input_fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
            current_sequence_header = ""
            current_sequence_lines = []
            
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    # New sequence header
                    if current_sequence_header:
                        # Process the previous sequence
                        exclude = False
                        for exclude_id in ids_to_exclude:
                            if exclude_id in current_sequence_header:
                                exclude = True
                                break
                        if not exclude:
                            outfile.write(current_sequence_header + '\n')
                            outfile.write('\n'.join(current_sequence_lines) + '\n')
                    
                    current_sequence_header = line
                    current_sequence_lines = []
                else:
                    # Sequence line
                    current_sequence_lines.append(line)
            
            # Process the last sequence
            if current_sequence_header:
                exclude = False
                for exclude_id in ids_to_exclude:
                    if exclude_id in current_sequence_header:
                        exclude = True
                        break
                if not exclude:
                    outfile.write(current_sequence_header + '\n')
                    outfile.write('\n'.join(current_sequence_lines) + '\n')

        print(f"Filtered sequences saved to '{output_file}'")

    except FileNotFoundError as e:
        print(f"Error: {e}. Make sure the input files exist.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    # Assuming input files are in the same directory as the script
    # and the output file will be created there.
    input_fasta = "clustered_core_proteins.faa"
    ids_csv = "human_hits.csv"
    output_txt = "filtered_sequences.txt"
    
    filter_sequences(input_fasta, ids_csv, output_txt)
