#!/usr/bin/env python3
"""Check the Python environment used by FAST XGB training."""

import importlib
import importlib.util
import sys
import warnings


DEPENDENCIES = (
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("scikit-learn", "sklearn"),
    ("xgboost", "xgboost"),
    ("optuna", "optuna"),
)
EXIT_MISSING = 10
EXIT_BROKEN = 11
EXIT_MISSING_AND_BROKEN = 12


def _one_line(value):
    return " ".join(str(value).replace("\t", " ").splitlines())


def check_dependencies():
    """Import every dependency independently and return results plus modules."""
    results = []
    modules = {}

    for package_name, import_name in DEPENDENCIES:
        try:
            available = importlib.util.find_spec(import_name) is not None
        except Exception as exc:  # pragma: no cover - unusual import hook failure
            results.append(
                {
                    "package": package_name,
                    "import_name": import_name,
                    "status": "broken",
                    "detail": "{}: {}".format(type(exc).__name__, _one_line(exc)),
                }
            )
            continue

        if not available:
            results.append(
                {
                    "package": package_name,
                    "import_name": import_name,
                    "status": "missing",
                    "detail": "",
                }
            )
            continue

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                module = importlib.import_module(import_name)
            modules[import_name] = module
            results.append(
                {
                    "package": package_name,
                    "import_name": import_name,
                    "status": "ok",
                    "detail": _one_line(getattr(module, "__version__", "unknown")),
                }
            )
        except Exception as exc:
            results.append(
                {
                    "package": package_name,
                    "import_name": import_name,
                    "status": "broken",
                    "detail": "{}: {}".format(type(exc).__name__, _one_line(exc)),
                }
            )

    return results, modules


def require_dependencies():
    """Return imported modules, or raise one error describing every problem."""
    results, modules = check_dependencies()
    missing = [result["package"] for result in results if result["status"] == "missing"]
    broken = [result for result in results if result["status"] == "broken"]

    problems = []
    if missing:
        problems.append("missing package(s): " + ", ".join(missing))
    if broken:
        problems.append(
            "package import failure(s): "
            + "; ".join(
                "{} ({})".format(result["package"], result["detail"])
                for result in broken
            )
        )
    if problems:
        raise RuntimeError("Python dependency check failed: " + ". ".join(problems) + ".")

    return modules


def print_report(results):
    missing = [result["package"] for result in results if result["status"] == "missing"]
    broken = [result for result in results if result["status"] == "broken"]

    if not missing and not broken:
        print("Python is ready for FAST XGBoost.")
        print("Using: {}".format(_one_line(sys.executable)))
        return

    print("Python is not ready for FAST XGBoost.")
    print("Using: {}".format(_one_line(sys.executable)))
    if missing:
        print("Missing required Python packages: {}.".format(", ".join(missing)))
    if broken:
        print(
            "Could not load these Python packages: {}.".format(
                ", ".join(result["package"] for result in broken)
            )
        )


def main():
    results, _ = check_dependencies()
    print_report(results)
    has_missing = any(result["status"] == "missing" for result in results)
    has_broken = any(result["status"] == "broken" for result in results)
    if has_missing and has_broken:
        return EXIT_MISSING_AND_BROKEN
    if has_missing:
        return EXIT_MISSING
    if has_broken:
        return EXIT_BROKEN
    return 0


if __name__ == "__main__":
    sys.exit(main())
