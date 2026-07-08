import os
import sys

import numpy as np

# Add the build lib directory to PYTHONPATH to find marmot extension
sys.path.insert(0, os.path.abspath("../../build/lib"))

import marmot

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
step.strainIncrementTarget = np.array([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
step.isStrainComponentControlled = np.array([True, True, True, True, True, True])
step.isStressComponentControlled = np.logical_not(step.isStrainComponentControlled)
step.stressIncrementTarget = np.zeros(6, dtype=np.float64)

solver.addStep(step)
solver.solve()

history = solver.getHistory()
if len(history) > 0:
    print(f"Final Stress 11: {history[-1].stress[0]}")
