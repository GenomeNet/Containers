#!/usr/bin/env python3

import argparse
import os
import subprocess
from pathlib import Path
import shutil
import tempfile
import pty
from termcolor import colored
import datetime

def is_fasta(filename):
    with open(filename, 'r') as file:
        first_line = file.readline().strip()
    return first_line.startswith(">")

def get_default_output_dir():
    date_str = datetime.date.today().strftime("deepG_virus_%Y%m%d_1")
    output_path = Path(date_str)
    counter = 1
    while output_path.exists():
        counter += 1
        date_str = datetime.date.today().strftime(f"deepG_virus_%Y%m%d_{counter}")
        output_path = Path(date_str)
    return date_str

def is_podman_installed():
    try:
        subprocess.run(["podman", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        return False
        
def run_podman(input_file, output_dir, use_gpu=False, batch_size=3000, level="binary", verbose=False, by_entry=False, fast=False):
    if not is_podman_installed():
        print("Error: Podman is not installed.")
        print("Please follow the installation instructions at https://podman.io/getting-started/installation.html")
        return

    input_path = Path(input_file).resolve()
    output_path = Path(output_dir).resolve()

    # Check if the input file is in the FASTA format
    if not is_fasta(input_path):
        print("Error: Input file is not in the FASTA format.")
        return

    with tempfile.TemporaryDirectory() as temp_data_dir:
        temp_data_path = Path(temp_data_dir)

        # Copy input file to temporary data directory
        shutil.copy(input_path, temp_data_path / input_path.name)

        cmd = [
            "podman", "run", "--rm", "-i",
            "--runtime", "nvidia" if use_gpu else "runc",
            "-v", f"{temp_data_path}:/data",
            "genomenet/virus",
            "predict.r",
            "--model", level,
            "--batch_size", str(batch_size),
            "--input", f"/data/{input_path.name}",
        ]
        
        if by_entry:
            cmd.append("--by_entry")
        if fast:
            cmd.append("--fast")

        print(colored("Starting deepG container", 'blue'))

        def read_output(fd):
            while True:
                try:
                    data = os.read(fd, 1024)
                except OSError:
                    break

                if not data:
                    break

                for line in data.decode().splitlines():
                    print(f"deepG: {line}")

        # Create a pseudo-terminal and run the command
        pid, fd = pty.fork()

        if pid == 0:
            os.execvp(cmd[0], cmd)
        else:
            read_output(fd)
            _, exit_status = os.waitpid(pid, 0)

        if exit_status != 0:
            print("Error: Container command failed.")
            success = False
        else:
            success = True

        if success:
            # Create the output directory if it doesn't exist
            if not output_path.exists():
                output_path.mkdir(parents=True)

        # Copy output files from temporary data directory to the specified output directory
        output_files = (temp_data_path / "output").glob("*")
        for output_file in output_files:
            if output_file.is_file():
                if verbose:
                    print(f"Copying output file '{output_file.name}' to the specified output directory...")
                shutil.copy(output_file, output_path / output_file.name)
            elif output_file.is_dir():  # Add this block to handle directories
                if verbose:
                    print(f"Copying output directory '{output_file.name}' to the specified output directory...")
                shutil.copytree(output_file, output_path / output_file.name)

        # Print a finishing message with input and output files information
        finishing_message = f"deepG container finished successfully.\nInput: '{input_path}', Output directory: '{output_path}'"
        print(colored(finishing_message, 'green'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Container with input and output paths",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example usage:\n"
               "  ./virus.py --input /path/to/example.fasta --output /path/to/output_dir\n"
               "  ./virus.py --input /path/to/example.fasta --output /path/to/output_dir --gpu\n"
               "  ./virus.py --input /path/to/example.fasta --output /path/to/output_dir --verbose"
    )
    parser.add_argument("--input", required=True, help="Path to input file")
    parser.add_argument("--output", default=get_default_output_dir(), help="Path to output directory (default: current date, e.g., 'deepG_virus_20210312_1')")
    parser.add_argument("--by_entry", action="store_true", help="Get predictions for each FASTA entry")
    parser.add_argument("--fast", action="store_true", help="Only evaluate one subsample per entry")
    parser.add_argument("--gpu", action="store_true", help="Use GPU mode")
    parser.add_argument("--batch_size", type=int, default=3000, help="Batch size for processing (default: 3000)")
    parser.add_argument("--level", choices=["binary", "genus"], default="binary", help="Prediction level (default: 'binary')")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose mode for more information during execution")

    args = parser.parse_args()
    run_podman(args.input, args.output, args.gpu, args.batch_size, args.level, args.verbose, args.by_entry, args.fast)