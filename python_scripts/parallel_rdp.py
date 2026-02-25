#!/usr/bin/env python3
# Alex Song, July 2025
# Used to run the rdp classifer in parallel
import argparse
import concurrent.futures
import math
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO


def split_fasta(input_file, num_chunks, output_prefix):
    """
    Splits the input FASTA file into num_chunks parts using FASTA records.
    Returns a list of chunk file paths.
    """

    # Parse the input FASTA file
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} not found")

    try:
        records = []
        with open(input_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                records.append(record)

        if not records:
            raise ValueError(f"No valid FASTA records found in {input_file}")

        total_records = len(records)
        chunk_size = math.ceil(total_records / num_chunks)
    except Exception as e:
        raise RuntimeError(f"Error parsing FASTA file {input_file}: {str(e)}")
    chunk_files = []

    for i in range(num_chunks):
        # Select records for this chunk
        chunk_records = records[i * chunk_size : (i + 1) * chunk_size]
        chunk_filename = f"{output_prefix}_chunk_{i}.fasta"
        with open(chunk_filename, "w") as f:
            SeqIO.write(chunk_records, f, "fasta")
        chunk_files.append(chunk_filename)
    return chunk_files


def resolve_classifier_path(classifier_path):
    """
    Resolve the classifier properties file path, allowing common relative locations.
    Returns an absolute path if found, otherwise raises FileNotFoundError.
    """
    path = Path(classifier_path).expanduser()
    candidates = []

    if path.is_absolute():
        candidates.append(path)
    else:
        cwd = Path.cwd()
        candidates.extend(
            [
                cwd / path,
                cwd / "runtime" / "classifiers" / path,
                cwd / "config" / "classifiers" / path,
            ]
        )

    checked = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_file():
            return str(resolved)
        checked.append(str(resolved))

    searched = "\n  - ".join(checked) if checked else str(path)
    raise FileNotFoundError(
        "RDP classifier properties file not found for '-t'.\n"
        f"Requested path: {classifier_path}\n"
        f"Checked paths:\n  - {searched}\n"
        "Provide an absolute path or place the classifier under runtime/classifiers/."
    )


def normalize_options(options):
    """
    Parse and normalize RDP options once before parallel execution.
    If '-t' is provided, validate and rewrite it to an absolute file path.
    """
    option_tokens = shlex.split(options)
    if "-t" in option_tokens:
        idx = option_tokens.index("-t")
        if idx + 1 >= len(option_tokens):
            raise ValueError("Invalid --options: '-t' requires a classifier properties file path")
        option_tokens[idx + 1] = resolve_classifier_path(option_tokens[idx + 1])
    return option_tokens


def run_rdp_classifier(chunk_file, memory_flag, option_tokens, result_file, timeout=43200):
    """
    Executes the RDP Classifier on a given chunk file.
    Command format:
        rdp_classifier <memory_flag> classify <options> -o <result_file> <chunk_file>
    A timeout (in seconds) is applied to prevent hanging.
    """
    cmd = (
        ["rdp_classifier", f"-Xmx{memory_flag}", "classify"]
        + option_tokens
        + ["-o", result_file, chunk_file]
    )
    print(f"Running command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
        if result.returncode != 0:
            print(f"Error processing {chunk_file}:\n{result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
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
    parser = argparse.ArgumentParser(
        description="Parallelize RDP Classifier execution on a FASTA file."
    )
    parser.add_argument("--input", required=True, help="Path to the input FASTA file")
    parser.add_argument(
        "--output", required=True, help="Path for the final output file"
    )
    parser.add_argument(
        "--threads", type=int, default=4, help="Number of parallel threads to use"
    )
    parser.add_argument(
        "--memory", required=True, help="Memory flag for rdp_classifier (e.g., '10g')"
    )
    parser.add_argument(
        "--options",
        required=True,
        help="Additional options for rdp_classifier (e.g., '-t /path/to/rRNAClassifier.properties')",
    )
    args = parser.parse_args()
    option_tokens = normalize_options(args.options)

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
            executor.submit(
                run_rdp_classifier, chunk, args.memory, option_tokens, res_file
            ): chunk
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
