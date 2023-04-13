import shutil
from argparse import ArgumentParser
from common import update_overview
from pathlib import Path

def main() -> None:
    parser = ArgumentParser(description="Duplicate a solver.")
    parser.add_argument("source", type=str, help="the name of the solver to duplicate")
    parser.add_argument("target", type=str, help="the name of the new solver")

    args = parser.parse_args()

    source_cpp_file = Path(__file__).parent.parent / "src" / f"{args.source}.cpp"
    if not source_cpp_file.is_file():
        raise RuntimeError(f"{source_cpp_file} does not exist")

    target_cpp_file = Path(__file__).parent.parent / "src" / f"{args.target}.cpp"
    if target_cpp_file.is_file():
        raise RuntimeError(f"{target_cpp_file} already exists")

    source_results_directory = Path(__file__).parent.parent / "results" / "output" / args.source
    if source_results_directory.is_dir():
        target_results_directory = Path(__file__).parent.parent / "results" / "output" / args.target
        if target_results_directory.is_dir():
            raise RuntimeError(f"{target_results_directory} already exists")

        print(f"Copying {source_results_directory} to {target_results_directory}")
        shutil.copytree(source_results_directory, target_results_directory)

    print(f"Copying {source_cpp_file} to {target_cpp_file}")
    shutil.copy(source_cpp_file, target_cpp_file)

    update_overview()

if __name__ == "__main__":
    main()
