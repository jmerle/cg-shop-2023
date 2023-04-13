import zipfile
from pathlib import Path

def main() -> None:
    results_directory = Path(__file__).parent.parent / "results"

    solution_files = []
    for instance_file in (results_directory / "instances").glob("*.instance.json"):
        instance_name = instance_file.name.split(".instance.json")[0]

        score_files = (results_directory / "output").rglob(f"{instance_name}.txt")
        best_score_file = min(score_files, key=lambda file: int(file.read_text(encoding="utf-8")))

        solution_files.append(best_score_file.parent / f"{instance_name}.json")

    submission_file = results_directory / "submission.zip"
    if submission_file.is_file():
        submission_file.unlink()

    with zipfile.ZipFile(submission_file, "w", zipfile.ZIP_DEFLATED) as file:
        for solution_file in solution_files:
            file.write(solution_file, solution_file.name)

    print(f"Successfully created submission archive at {submission_file}")

if __name__ == "__main__":
    main()
