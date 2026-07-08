import os
import numpy as np

def compareHistoryToReferenceSolution(history, reference_file):
    """
    Compare a solver history against a reference npz file.
    If the reference file does not exist, it creates one.
    """
    if len(history) == 0:
        raise RuntimeError("Solver history is empty!")

    stress = np.array([h.stress for h in history])
    if hasattr(history[0], 'F'):
        strain = np.array([h.F for h in history])
    else:
        strain = np.array([h.strain for h in history])

    if not os.path.exists(reference_file):
        np.savez(reference_file, stress=stress, strain=strain)
        print(f"Reference solution created: {reference_file}")
    else:
        ref = np.load(reference_file)
        np.testing.assert_allclose(stress, ref['stress'], atol=1e-6, rtol=1e-5)
        np.testing.assert_allclose(strain, ref['strain'], atol=1e-6, rtol=1e-5)
        print(f"Reference solution passed: {reference_file}")
