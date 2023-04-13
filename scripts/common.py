import json
from pathlib import Path

def update_overview() -> None:
    scores_by_solver = {}
    outputs_root = Path(__file__).parent.parent / "results" / "output"

    for directory in outputs_root.iterdir():
        scores_by_instance = {}

        for file in directory.iterdir():
            if file.name.endswith(".txt"):
                scores_by_instance[file.stem] = float(file.read_text(encoding="utf-8").strip())

        scores_by_solver[directory.name] = scores_by_instance

    overview_template_file = Path(__file__).parent.parent / "results" / "overview.tmpl.html"
    overview_file = Path(__file__).parent.parent / "results" / "overview.html"

    overview_template = overview_template_file.read_text(encoding="utf-8")
    overview = overview_template.replace("/* scores_by_solver */{}", json.dumps(scores_by_solver))

    with overview_file.open("w+", encoding="utf-8") as file:
        file.write(overview)

    print(f"Overview: file://{overview_file.resolve()}")
