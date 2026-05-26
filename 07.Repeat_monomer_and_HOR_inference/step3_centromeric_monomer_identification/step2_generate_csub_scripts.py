--- START OF FILE 02.generate_csub_scripts.py ---

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import math

def generate_scripts(sh_file, csub_header_info, chunk_size, output_dir):
    """
    Main function to generate csub job submission scripts by splitting a task list.
    """
    # 1. Parse csub header parameters
    print("--- Step 1/4: Parsing csub header parameters ---")
    header_parts = csub_header_info.split(';')
    if len(header_parts) < 3:
        print(f"Error: Argument 2 (csub header) is incorrectly formatted. At least 3 values separated by ';' are required.", file=sys.stderr)
        print("Example: 'job_name;8;/path/to/workdir'", file=sys.stderr)
        sys.exit(1)
        
    job_name = header_parts[0]
    num_cores = header_parts[1]
    work_dir = header_parts[2]
    
    print(f"  - Job Name (-J, -o, -e): {job_name}")
    print(f"  - Core Count (-n): {num_cores}")
    print(f"  - Working Directory (-cwd, cd): {work_dir}")

    # Define the csub header template
    # Note: Including %J (Job ID) in output/error files prevents overwriting across multiple runs.
    csub_template = f"""#!/bin/bash
#CSUB -J {job_name}
#CSUB -q c01
#CSUB -o {job_name}.%J.o
#CSUB -e {job_name}.%J.e
#CSUB -n {num_cores}
#CSUB -R span[hosts=1]
#CSUB -cwd {work_dir}
cd {work_dir}

"""

    # 2. Read the list of tasks to be analyzed
    print("\n--- Step 2/4: Reading task file ---")
    try:
        with open(sh_file, 'r') as f:
            # Filter out empty lines
            tasks = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Task file '{sh_file}' not found.", file=sys.stderr)
        sys.exit(1)

    if not tasks:
        print(f"Warning: Task file '{sh_file}' is empty. No scripts were generated.", file=sys.stderr)
        sys.exit(0)
    
    print(f"Successfully loaded {len(tasks)} tasks.")

    # 3. Prepare the output directory
    print(f"\n--- Step 3/4: Preparing output directory '{output_dir}' ---")
    try:
        os.makedirs(output_dir, exist_ok=True)
        print("Output directory is ready.")
    except OSError as e:
        print(f"Error: Unable to create output directory '{output_dir}': {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Split tasks and generate scripts
    print(f"\n--- Step 4/4: Generating scripts with {chunk_size} tasks per file ---")
    
    num_chunks = math.ceil(len(tasks) / chunk_size)
    
    for i in range(num_chunks):
        start_index = i * chunk_size
        end_index = start_index + chunk_size
        task_chunk = tasks[start_index:end_index]
        
        # Define output filename
        output_filename = f"submit_part_{i+1}.sh"
        output_filepath = os.path.join(output_dir, output_filename)
        
        try:
            with open(output_filepath, 'w') as f_out:
                # Write header
                f_out.write(csub_template)
                # Write the current batch of tasks
                for task in task_chunk:
                    f_out.write(task + '\n')
            
            print(f"  - Generated: {output_filepath} (contains {len(task_chunk)} tasks)")
        
        except IOError as e:
            print(f"Error: Failed to write to file '{output_filepath}': {e}", file=sys.stderr)
            # Stop execution on IO error to prevent incomplete results
            sys.exit(1)
            
    print(f"\n--- Completion ---")
    print(f"Successfully generated {num_chunks} csub scripts in directory '{output_dir}'.")
    print("Batch submission command: 'for f in $(ls *.sh); do csub $f; done'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a shell file containing multiple commands into several submission scripts with csub headers.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "sh_file",
        help="Input file containing commands to execute (one task per line)."
    )
    parser.add_argument(
        "csub_header_info",
        help="csub header parameters, semicolon-separated (';').\n"
             "Format: 'job_name;cores;work_dir'\n"
             "Example: 'moddotplot;8;/share/project/analysis'"
    )
    parser.add_argument(
        "chunk_size",
        type=int,
        help="Number of tasks to include in each generated script."
    )
    parser.add_argument(
        "output_dir",
        help="Directory where the generated csub scripts will be stored."
    )

    args = parser.parse_args()

    generate_scripts(args.sh_file, args.csub_header_info, args.chunk_size, args.output_dir)