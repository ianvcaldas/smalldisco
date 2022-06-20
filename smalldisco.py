import argparse
import subprocess
import sys
from pathlib import Path


def get_command_line_arguments():
    parser = argparse.ArgumentParser("Locate siRNA regions and get small RNA tails.")
    parser.add_argument(
        "mode", type=str, help="Analysis mode. Can be 'align', 'locate', 'tail', or 'run_all'."
    )
    return parser.parse_known_args()


def run_snakemake(my_args, extra_args):
    targets = {
        "align": "align_all_samples",
        "sirna": "locate_sirna_regions",
        "tail": "get_tailing",
        "run_all": "get_tailing"
    }
    try:
      target = targets[my_args.mode]
    except KeyError:
        raise Exception(f"Unrecognized mode: {my_args.mode}")
    command = [
        "snakemake",
        "-c1",
        "--snakefile",
        "./workflow/Snakefile"
    ]
    command.extend(extra_args)
    command.append("--force")
    command.append(target)
    print(" ".join(command))
    subprocess.run(command)


if __name__ == "__main__":
    my_args, extra_args = get_command_line_arguments()
    run_snakemake(my_args, extra_args)
