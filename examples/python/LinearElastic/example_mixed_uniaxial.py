import marmot
import numpy as np

print("Running example for HypoElastic material: LINEARELASTIC")

# Material properties extracted from C++ tests
properties = np.array([20000, 0.25], dtype=np.float64)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("LINEARELASTIC", properties, options)

# Setup a loading step
step = marmot.solvers.HypoElasticSolver.Step()
step.timeStart = 0.0
step.timeEnd = 1.0
step.dTStart = 0.1

step.isStrainComponentControlled = np.array([True, False, False, True, True, True], dtype=bool)
step.isStressComponentControlled = np.array([False, True, True, False, False, False], dtype=bool)

strain_target = np.zeros(6, dtype=np.float64)
strain_target[0] = 0.1
step.strainIncrementTarget = strain_target

stress_target = np.zeros(6, dtype=np.float64)
step.stressIncrementTarget = stress_target

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
