import glob
import os
import argparse
import gzip

# usage "python convert_ml_fasta.py "*fasta"

def convert_multiline_fasta_to_single_line(input_fasta, output_fasta):
    # Determine if the file is gzipped
    is_gzipped = input_fasta.endswith(".gz")
    
    # Open the input file with appropriate method
    open_func = gzip.open if is_gzipped else open
    
    with open_func(input_fasta, 'rt') as infile, open(output_fasta, 'w') as outfile:
        sequence = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    outfile.write(sequence + "\n")
                    sequence = ""
                outfile.write(line + "\n")
            else:
                sequence += line
        if sequence:
            outfile.write(sequence + "\n")

def process_files(pattern):
    # Use glob to find all files matching the pattern
    fasta_files = glob.glob(pattern)
    
    for input_fasta in fasta_files:
        # Define the output file name
        output_fasta = f"converted_{os.path.basename(input_fasta).replace('.gz', '')}"
        print(f"Processing {input_fasta} -> {output_fasta}")
        # Convert the file
        convert_multiline_fasta_to_single_line(input_fasta, output_fasta)

def main():
    parser = argparse.ArgumentParser(description="Convert multiline FASTA files to single-line format.")
    parser.add_argument("pattern", type=str, help="Wildcard pattern to match input FASTA files.")
    args = parser.parse_args()

    process_files(args.pattern)

if __name__ == "__main__":
    main()