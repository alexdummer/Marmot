import argparse
import os

import numpy as np


def writeReferenceSolution(history, reference_file):
    """Write solver history to a reference npz file."""
    if len(history) == 0:
        raise RuntimeError("Solver history is empty!")

    stress = np.array([h.stress for h in history])
    if hasattr(history[0], "F"):
        strain = np.array([h.F for h in history])
    else:
        strain = np.array([h.strain for h in history])

    np.savez(reference_file, stress=stress, strain=strain)
    print(f"Reference solution created: {reference_file}")


def compareHistoryToReferenceSolution(history, reference_file):
    """Compare a solver history against a reference npz file."""
    if len(history) == 0:
        raise RuntimeError("Solver history is empty!")

    if not os.path.exists(reference_file):
        raise FileNotFoundError(
            f"Reference file not found: {reference_file}. Run with --writeReferenceSolution to create it."
        )

    stress = np.array([h.stress for h in history])
    if hasattr(history[0], "F"):
        strain = np.array([h.F for h in history])
    else:
        strain = np.array([h.strain for h in history])

    ref = np.load(reference_file)
    np.testing.assert_allclose(stress, ref["stress"], atol=1e-6, rtol=1e-5)
    np.testing.assert_allclose(strain, ref["strain"], atol=1e-6, rtol=1e-5)
    print(f"Reference solution passed: {reference_file}")


def run_test(history, script_file):
    """
    Unified entry point for testing python scripts.
    Parses command line arguments to either compare or write the reference solution.
    """
    parser = argparse.ArgumentParser(description="Marmot example test script")
    parser.add_argument("--writeReferenceSolution", action="store_true", help="Write the reference solution")

    args, _ = parser.parse_known_args()

    reference_file = os.path.splitext(script_file)[0] + "_reference.npz"

    if args.writeReferenceSolution:
        writeReferenceSolution(history, reference_file)
    else:
        compareHistoryToReferenceSolution(history, reference_file)
