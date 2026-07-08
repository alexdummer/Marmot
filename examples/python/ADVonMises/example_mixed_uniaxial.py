import marmot
import numpy as np

print("Running example for HypoElastic material: ADVONMISES")

# Material properties extracted from C++ tests
properties = np.array([210000.0, 0.3, 200.0, 2100.0, 20.0, 20], dtype=np.float64)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("ADVONMISES", properties, options)

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
if len(history) > 0:
    print(f"Final Stress 11: {history[-1].stress[0]}\nFinal Strain 11: {history[-1].strain[0]}\nFinal Strain 22: {history[-1].strain[1]}")
