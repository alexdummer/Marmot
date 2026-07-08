import os
import sys

import numpy as np

# Add the build lib directory to PYTHONPATH to find marmot extension
sys.path.insert(0, os.path.abspath("../../build/lib"))

import marmot

print("Running example for HypoElastic material: LINEARVISCOELASTICPOWERLAW")

# Material properties extracted from C++ tests
properties = np.array([2e5, 0.2, 0.5, 0.1, 10.0, 0.0001, 1.0], dtype=np.float64)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("LINEARVISCOELASTICPOWERLAW", properties, options)

# Setup a loading step
step = marmot.solvers.HypoElasticSolver.Step()
step.strainIncrementTarget = np.array([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
step.isStrainComponentControlled = np.array([True, True, True, True, True, True])
step.isStressComponentControlled = np.logical_not(step.isStrainComponentControlled)
step.stressIncrementTarget = np.zeros(6, dtype=np.float64)

solver.addStep(step)
solver.solve()

history = solver.getHistory()
if len(history) > 0:
    print(f"Final Stress 11: {history[-1].stress[0]}")
