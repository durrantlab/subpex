"""Simple wrapper around bash script."""

import argparse
import os
import subprocess

current_dir = os.path.dirname(os.path.realpath(__file__))

# Path to the Bash script relative to the current Python script
bash_script_path = os.path.join(current_dir, "sp_restart.sh")


def main():
    parser = argparse.ArgumentParser(
        description="Restart SubPEx simulation with options."
    )
    parser.add_argument(
        "-n", "--generation", type=int, help="Generation number for truncation"
    )
    parser.add_argument("-e", "--env", type=str, help="Conda environment to activate")
    args = parser.parse_args()

    # Build the command to run the Bash script with the necessary arguments
    command = [bash_script_path]
    if args.generation is not None:
        command.extend(["-n", str(args.generation)])
    if args.env is not None:
        command.extend(["-e", str(args.env)])

    # Execute the Bash script
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
