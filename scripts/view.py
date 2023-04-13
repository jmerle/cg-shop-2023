import argparse
import matplotlib.pyplot as plt
import os
import subprocess
import tempfile
from adjustText import adjust_text
from cgshop2023_pyutils import read_instance, read_solution
from contextlib import contextmanager
from matplotlib.collections import PolyCollection
from pathlib import Path
from typing import Generator

@contextmanager
def temporary_file(suffix: str) -> Generator[Path, None, None]:
    fd, name = tempfile.mkstemp(suffix=suffix)
    os.close(fd)

    path = Path(name)

    try:
        yield path
    finally:
        path.unlink()

def main() -> None:
    parser = argparse.ArgumentParser(description="View a solver's solution for an instance.")
    parser.add_argument("solver", type=str, help="the solver to show the solution of")
    parser.add_argument("instance", type=str, help="the instance to show the solution of")
    parser.add_argument("--coordinates", action="store_true", help="show coordinates of points in instance polygon")

    args = parser.parse_args()

    instance_file = Path(__file__).parent.parent / "results" / "instances" / f"{args.instance}.instance.json"
    if not instance_file.is_file():
        raise RuntimeError(f"Instance not found, {instance_file} is not a file")

    solution_file = Path(__file__).parent.parent / "results" / "output" / args.solver / f"{args.instance}.json"
    if not solution_file.is_file():
        raise RuntimeError(f"Solution not found, {solution_file} is not a file")

    instance = read_instance(instance_file)
    solution = read_solution(solution_file)

    fig = plt.figure(figsize=(50, 50))
    ax = plt.gca()

    if args.coordinates:
        texts = []
        for polygon in [instance["outer_boundary"]] + instance["holes"]:
            for point in polygon:
                x = point["x"]
                y = point["y"]
                texts.append(ax.text(x, y, f"({x}, {y})"))

    for parsed_polygons, color, style in [
        ([instance["outer_boundary"]] + instance["holes"], "#000000", "-"),
        (solution["polygons"], "#FF0000", ":")
    ]:
        polygons = []
        for polygon in parsed_polygons:
            polygons.append([[point["x"], point["y"]] for point in polygon])

        ax.add_collection(PolyCollection(polygons, edgecolors=color, facecolors="None", linestyles=style))

    ax.autoscale_view()

    if args.coordinates:
        adjust_text(texts)

    with temporary_file(".png") as file:
        fig.savefig(file, bbox_inches="tight", pad_inches=0)
        subprocess.run(["qview", str(file)])

if __name__ == "__main__":
    main()
