"""The upload workflow itself gates the release, without running tests twice.

Static contract: follow required jobs and local scripts before the upload.
This is not a proof of a live GitHub response; the exact implementation and
failure handling must also be reviewed. No API calls or publishing occur here.
"""

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/publish.yml"


def _workflow():
    return yaml.safe_load(WORKFLOW.read_text())


def _gate_steps():
    jobs = _workflow()["jobs"]
    uploads = [
        (name, i)
        for name, job in jobs.items()
        for i, step in enumerate(job.get("steps", []))
        if "pypi-publish@" in step.get("uses", "")
    ]
    assert len(uploads) == 1, "One identifiable PyPI upload must be guarded"
    job_name, index = uploads[0]
    seen = set()
    steps = []

    def ancestors(name):
        if name in seen:
            return
        seen.add(name)
        job = jobs[name]
        assert not job.get("continue-on-error"), "A failed gate must block publishing"
        assert "always()" not in str(job.get("if", "")), "Do not bypass failed required jobs"
        needs = job.get("needs", [])
        for dependency in [needs] if isinstance(needs, str) else needs:
            ancestors(dependency)
        steps.extend(job.get("steps", []) if name != job_name else job["steps"][:index])

    ancestors(job_name)
    upload = jobs[job_name]["steps"][index]
    assert "always()" not in str(upload.get("if", ""))
    return steps


def _gate_text():
    steps = _gate_steps()
    text = yaml.safe_dump(steps)
    # Permit an inline check or a checked-in helper actually called by a step.
    # Resolve references transitively, never scanning unrelated repo scripts.
    seen = set()
    while True:
        paths = set(re.findall(r"(?:\./)?(?:scripts|\.github)/[\w./-]+\.(?:py|sh|js)", text))
        paths -= seen
        if not paths:
            break
        for name in paths:
            seen.add(name)
            path = ROOT / name
            assert path.is_file(), f"Missing release-gate helper: {name}"
            text += "\n" + path.read_text()
    return text


def test_publish_checks_out_the_release_tag_not_the_default_branch():
    checkouts = [s for s in _gate_steps() if s.get("uses", "").startswith("actions/checkout@")]
    assert checkouts, "The checked sources must be the release's tagged sources"
    assert any(
        "github.event.release.tag_name" in str(s.get("with", {}).get("ref", "")) for s in checkouts
    ), "Explicitly check out the release tag being published"


def test_publish_compares_tag_with_pyproject_package_version():
    text = _gate_text()
    assert "pyproject.toml" in text, "Read the version of the package being built"
    assert re.search(r"(?:tag_name|GITHUB_REF|release.tag|RELEASE_TAG)", text)
    assert "version" in text and re.search(r"(?:tomllib|toml|project)", text)
    assert re.search(r"(?:!=|==| -ne | -eq | -n?z |assert )", text)
    assert re.search(r"(?:exit [1-9]|raise |setFailed|sys.exit)", text), (
        "Mismatch must stop upload"
    )


def test_publish_requires_the_ci_workflow_not_any_successful_workflow():
    assert "ci.yml" in _gate_text(), "Require the actual CI workflow, not an unrelated green run"


def test_publish_selects_ci_for_the_tagged_commit():
    text = _gate_text()
    assert re.search(r"(?:head_sha|--commit\b)", text), "CI must belong to the tagged commit"
    assert re.search(r"(?:git rev-parse|GITHUB_SHA|github.sha|target_commitish)", text)


def test_publish_requires_completed_success_and_fails_closed():
    text = _gate_text()
    assert "success" in text and "conclusion" in text, "A queued or failed CI run is not green"
    assert re.search(r"(?:completed|gh run watch|wait_for|waitFor)", text)
    assert re.search(r"(?:exit [1-9]|raise |setFailed|sys.exit)", text)
    assert "|| true" not in text, "A missing/error CI response must not authorize upload"


def test_publish_smokes_the_built_wheel_in_a_clean_venv():
    text = _gate_text()
    assert re.search(r"(?:-m venv|virtualenv)", text), "Install the wheel in a clean environment"
    assert "--system-site-packages" not in text
    assert re.search(r"(?:pip install|pip.*-m| -m pip)", text)
    assert re.search(r"(?:dist/.*whl|\.whl|wheel_path|wheel_file)", text), (
        "Install the built artifact"
    )
    assert not re.search(r"pip install[^\n]*(?: -e | --editable )", text)


def test_publish_checks_installed_wheel_version_before_upload():
    text = _gate_text()
    assert "--version" in text, "The installed wheel must report the release version"
    assert re.search(
        r"(?:bin/allelix|Scripts/allelix|activate|venv_python|venv.*-m allelix)", text
    )
    assert "version" in text and re.search(r"(?:!=|==|assert | -eq | -ne )", text)


def test_publish_does_not_run_the_test_suite_again():
    text = _gate_text()
    assert not re.search(r"(?:\bpytest\b|\btox\b|unittest discover)", text)
    assert not any("/ci.yml" in str(job.get("uses", "")) for job in _workflow()["jobs"].values())


def test_preupload_steps_cannot_ignore_gate_errors():
    for step in _gate_steps():
        assert not step.get("continue-on-error"), step
        assert "always()" not in str(step.get("if", "")), step
