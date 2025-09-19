import sys

def filter_fasta_by_pattern(fasta_file, pattern_file, output_file):
    """
    Filters sequences from a FASTA file whose headers partially match
    a list of patterns from a CSV file, and writes them to an output file.

    Args:
        fasta_file (str): Path to the input FASTA file.
        pattern_file (str): Path to the input pattern CSV file.
        output_file (str): Path to the output FASTA file where results will be saved.
    """
    # Read patterns from the pattern CSV file
    patterns = set()
    try:
        with open(pattern_file, 'r') as f_pattern:
            for line in f_pattern:
                pattern = line.strip()
                if pattern:
                    patterns.add(pattern)
    except FileNotFoundError:
        print(f"Error: Pattern file '{pattern_file}' not found.", file=sys.stderr)
        return

    if not patterns:
        print("Warning: No patterns found in the pattern file. No sequences will be filtered.", file=sys.stderr)
        return

    # Process the FASTA file and write to output
    current_header = None
    current_sequence_lines = []
    
    try:
        with open(fasta_file, 'r') as f_fasta, open(output_file, 'w') as f_out:
            for line in f_fasta:
                line = line.strip()
                if line.startswith('>'):
                    # If we have a sequence block, check for a match
                    if current_header:
                        for pattern in patterns:
                            if pattern in current_header:
                                f_out.write(current_header + '\n')
                                f_out.write("".join(current_sequence_lines) + '\n')
                                break # Move to the next header once a match is found for the current sequence
                    
                    # Start a new sequence block
                    current_header = line
                    current_sequence_lines = []
                else:
                    # Append sequence lines
                    current_sequence_lines.append(line)
            
            # After the loop, check the last sequence block
            if current_header and current_sequence_lines:
                for pattern in patterns:
                    if pattern in current_header:
                        f_out.write(current_header + '\n')
                        f_out.write("".join(current_sequence_lines) + '\n')
                        break
        print(f"Filtering complete. Matching sequences saved to '{output_file}'.")
                        
    except FileNotFoundError:
        print(f"Error: One or more files not found. Check '{fasta_file}' or '{output_file}'.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    # Define your file names
    fasta_input_file = "merged_proteins.faa"
    pattern_input_file = "pattern.csv"
    output_target_file = "core_proteins.faa" # New output file name
    
    # Call the function to perform the filtering and write results
    filter_fasta_by_pattern(fasta_input_file, pattern_input_file, output_target_file)
