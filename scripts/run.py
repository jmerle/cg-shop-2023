import argparse
import multiprocessing
import subprocess
import sys
from cgshop2023_pyutils import read_instance, read_solution, verify
from common import update_overview
from pathlib import Path
from tqdm.contrib.concurrent import process_map

def run_instance(options: tuple[Path, Path, Path]) -> int:
    solver, instance, output_directory = options
    instance_name = instance.name.split(".")[0]

    stdout_file = output_directory / f"{instance_name}.json"
    stderr_file = output_directory / f"{instance_name}.log"

    with instance.open("rb") as stdin, stdout_file.open("wb+") as stdout, stderr_file.open("wb+") as stderr:
            try:
                process = subprocess.run([str(solver)], stdin=stdin, stdout=stdout, stderr=stderr, timeout=60_000)

                if process.returncode != 0:
                    raise RuntimeError(f"Solver exited with status code {process.returncode} for instance {instance_name}")
            except subprocess.TimeoutExpired:
                raise RuntimeError(f"Solver timed out on instance {instance_name}")

    try:
        parsed_instance = read_instance(instance)
        parsed_solution = read_solution(stdout_file)

        error = verify(parsed_instance, parsed_solution)

        if error != "":
            raise ValueError(error)
    except Exception as e:
        raise RuntimeError(f"Solver provided invalid output for instance {instance_name}: {str(e)}")

    return len(parsed_solution["polygons"])

def run(solver: Path, instances: list[Path], output_directory: Path) -> None:
    if not output_directory.is_dir():
        output_directory.mkdir(parents=True)

    try:
        scores = process_map(run_instance,
                             [(solver, instance, output_directory) for instance in instances],
                             max_workers=max(1, multiprocessing.cpu_count() - 4))
    except Exception as e:
        print(f"\033[91m{str(e)}\033[0m")
        sys.exit(1)

    for i, instance in enumerate(instances):
        instance_name = instance.name.split(".")[0]

        print(f"{instance_name}: {scores[i]}")

        score_file = output_directory / f"{instance_name}.txt"
        with score_file.open("w+", encoding="utf-8") as file:
            file.write(str(scores[i]))

    if len(instances) > 1:
        print(f"Total score: {sum(scores)}")

def main() -> None:
    parser = argparse.ArgumentParser(description="Run a solver.")
    parser.add_argument("solver", type=str, help="the solver to run")
    parser.add_argument("--instance", type=str, help="the instance to run on (defaults to all instances)")

    args = parser.parse_args()

    solver = Path(__file__).parent.parent / "cmake-build-release" / "bin" / args.solver
    if not solver.is_file():
        raise RuntimeError(f"Solver not found, {solver} is not a file")

    output_directory = Path(__file__).parent.parent / "results" / "output" / args.solver

    if args.instance is not None:
        instances = [Path(__file__).parent.parent / "results" / "instances" / f"{args.instance}.instance.json"]
        if not instances[0].is_file():
            raise RuntimeError(f"Instance not found, {instances[0]} is not a file")
    else:
        instances = sorted((Path(__file__).parent.parent / "results" / "instances").glob("*.instance.json"))

    run(solver, instances, output_directory)
    update_overview()

if __name__ == "__main__":
    main()
