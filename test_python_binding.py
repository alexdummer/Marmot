import os
import sys

# Add the build lib directory to PYTHONPATH to find marmot extension
sys.path.insert(0, os.path.abspath("build/lib"))

import marmot
import numpy as np

print("Successfully imported marmot!")
print("Modules available:", dir(marmot))
print("Solvers available:", dir(marmot.solvers))

options = marmot.solvers.FiniteStrainSolver.SolverOptions()
props = np.array([210000.0, 0.3])

print("Creating FiniteStrainSolver...")
solver = marmot.solvers.FiniteStrainSolver("COMPRESSIBLENEOHOOKE", props, options)

step = marmot.solvers.FiniteStrainSolver.Step()
step.gradUIncrementTarget = np.array([[0.1, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

step.isGradUComponentControlled = np.array([[True, True, True], [True, False, True], [True, True, False]])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
step.stressIncrementTarget = np.zeros((3, 3))

solver.addStep(step)
print("Solving...")
solver.solve()

history = solver.getHistory()
print(f"Solved successfully. History size: {len(history)}")
if len(history) > 0:
    print(f"Final Stress XX: {history[-1].stress[0,0]}")
