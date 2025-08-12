#!/usr/bin/env python3
import argparse
import math
import os
import subprocess
import tempfile
import shutil
import concurrent.futures
from Bio import SeqIO

def split_fasta(input_file, num_chunks, output_prefix):
    """
    Splits the input FASTA file into num_chunks parts using FASTA records.
    Returns a list of chunk file paths.
    """
    # Parse the input FASTA file
    records = list(SeqIO.parse(input_file, "fasta"))
    total_records = len(records)
    chunk_size = math.ceil(total_records / num_chunks)
    chunk_files = []

    for i in range(num_chunks):
        # Select records for this chunk
        chunk_records = records[i * chunk_size : (i + 1) * chunk_size]
        chunk_filename = f"{output_prefix}_chunk_{i}.fasta"
        with open(chunk_filename, "w") as f:
            SeqIO.write(chunk_records, f, "fasta")
        chunk_files.append(chunk_filename)
    return chunk_files

def run_rdp_classifier(chunk_file, memory_flag, options, result_file, timeout=43200):
    """
    Executes the RDP Classifier on a given chunk file.
    Command format:
        rdp_classifier <memory_flag> classify <options> -o <result_file> <chunk_file>
    A timeout (in seconds) is applied to prevent hanging.
    """
    cmd = f"rdp_classifier -Xmx{memory_flag} classify {options} -o {result_file} {chunk_file}"
    print(f"Running command: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            print(f"Error processing {chunk_file}:\n{result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd)
    except subprocess.TimeoutExpired:
        print(f"Command timed out for {chunk_file}")
        raise
    return result_file

def concatenate_files(file_list, output_file):
    """Concatenates the contents of each file in file_list into output_file."""
    with open(output_file, "w") as outfile:
        for fname in file_list:
            with open(fname, "r") as infile:
                shutil.copyfileobj(infile, outfile)

def cleanup_files(file_list):
    """Removes all files in file_list."""
    for f in file_list:
        try:
            os.remove(f)
        except Exception as e:
            print(f"Error removing file {f}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Parallelize RDP Classifier execution on a FASTA file.")
    parser.add_argument("--input", required=True, help="Path to the input FASTA file")
    parser.add_argument("--output", required=True, help="Path for the final output file")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads to use")
    parser.add_argument("--memory", required=True, help="Memory flag for rdp_classifier (e.g., '10g')")
    parser.add_argument("--options", required=True, help="Additional options for rdp_classifier (e.g., '-t /path/to/rRNAClassifier.properties')")
    args = parser.parse_args()

    # Create a temporary directory to store chunk files
    temp_dir = tempfile.mkdtemp(prefix="rdp_chunks_")
    output_prefix = os.path.join(temp_dir, "chunk")

    print("Splitting input FASTA file into records...")
    chunk_files = split_fasta(args.input, args.threads, output_prefix)

    # Prepare result filenames corresponding to each chunk
    result_files = [f"{chunk_file}.result" for chunk_file in chunk_files]

    print("Running RDP Classifier in parallel on FASTA chunks...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_chunk = {
            executor.submit(run_rdp_classifier, chunk, args.memory, args.options, res_file): chunk
            for chunk, res_file in zip(chunk_files, result_files)
        }
        for future in concurrent.futures.as_completed(future_to_chunk):
            chunk = future_to_chunk[future]
            try:
                res = future.result()
                print(f"Finished processing {chunk}, result saved in {res}")
            except Exception as exc:
                print(f"Error processing chunk {chunk}: {exc}")
                raise exc

    print("Concatenating all result files...")
    concatenate_files(result_files, args.output)
    print(f"Final output written to {args.output}")

    print("Cleaning up temporary files...")
    cleanup_files(chunk_files + result_files)
    os.rmdir(temp_dir)
    print("Done.")

if __name__ == "__main__":
    main()
