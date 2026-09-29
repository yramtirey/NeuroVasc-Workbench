import json
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CASES_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "cases"
)


app = FastAPI(
    title="NeuroVasc Workbench API",
    version="0.1.0",
)


app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
def root():
    return {
        "name": "NeuroVasc Workbench API",
        "version": "0.1.0",
    }


@app.get("/api/cases")
def list_cases():

    if not CASES_DIR.exists():
        return []

    cases = []

    for case_dir in sorted(
        CASES_DIR.iterdir()
    ):

        report_path = (
            case_dir
            / "case_report.json"
        )

        if not report_path.exists():
            continue

        with open(report_path) as file:
            report = json.load(file)

        cases.append(
            {
                "case_id":
                    report["case_id"],

                "source_file":
                    report["source_file"],

                "vessel_count":
                    report["vessel_count"],
            }
        )

    return cases


@app.get("/api/cases/{case_id}")
def get_case(
    case_id: str,
):

    report_path = (
        CASES_DIR
        / case_id
        / "case_report.json"
    )

    if not report_path.exists():

        raise HTTPException(
            status_code=404,
            detail=(
                f"Case '{case_id}' not found."
            ),
        )

    with open(report_path) as file:
        return json.load(file)
@app.get(
    "/api/cases/{case_id}/geometry/{label}"
)
def get_vessel_geometry(
    case_id: str,
    label: int,
):

    geometry_path = (
        CASES_DIR
        / case_id
        / "geometry"
        / f"label_{label}.vtp"
    )

    if not geometry_path.exists():

        raise HTTPException(
            status_code=404,
            detail=(
                f"Geometry for label {label} "
                f"in case '{case_id}' not found."
            ),
        )

    return FileResponse(
        geometry_path,
        media_type="application/xml",
        filename=f"label_{label}.vtp",
    )


@app.get("/api/cases/{case_id}/profiles/{label}")
def get_vessel_profile(case_id: str, label: int):
    profile_path = CASES_DIR / case_id / "profiles" / f"label_{label}.json"
    if not profile_path.is_file():
        raise HTTPException(
            status_code=404,
            detail=f"Profile for label {label} in case '{case_id}' not found.",
        )
    with profile_path.open() as file:
        return json.load(file)
