import marmot
import numpy as np

print("Running example for HypoElastic material: B4")

# Material properties extracted from C++ tests
properties = np.array(
    [
        0.2,
        20.6,
        91.0,
        4.8,
        5.9,
        0.1,
        0.5,
        12.0,
        1e-5,
        0.0,
        3.0,
        1.45,
        -4.5,
        -0.0015,
        90.0,
        7.0,
        1.0,
        400.0,
        11.0,
        2e-4,
        -100.0,
        1.0,
    ],
    dtype=np.float64,
)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("B4", properties, options)

# Setup a loading step
step = marmot.solvers.HypoElasticSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1
step.strainIncrementTarget = np.zeros(6, dtype=np.float64)
step.isStrainComponentControlled = np.array([False, False, False, True, True, True])
step.isStressComponentControlled = np.logical_not(step.isStrainComponentControlled)
s = np.zeros(6, dtype=np.float64)
s[0] = 10.0
step.stressIncrementTarget = s

solver.addStep(step)
solver.solve()

history = solver.getHistory()
if len(history) == 0:
    raise RuntimeError("Solver history is empty!")

import os

reference_file = os.path.splitext(__file__)[0] + "_reference.npz"

stress = np.array([h.stress for h in history])
if hasattr(history[0], "F"):
    strain = np.array([h.F for h in history])
else:
    strain = np.array([h.strain for h in history])

if not os.path.exists(reference_file):
    np.savez(reference_file, stress=stress, strain=strain)
    print(f"Reference solution created: {reference_file}")
else:
    ref = np.load(reference_file)
    np.testing.assert_allclose(stress, ref["stress"], atol=1e-6, rtol=1e-5)
    np.testing.assert_allclose(strain, ref["strain"], atol=1e-6, rtol=1e-5)
    print(f"Reference solution passed: {reference_file}")
