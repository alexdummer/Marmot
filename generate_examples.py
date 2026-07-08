import os
import re

hypoelastic_materials = [
    "VONMISES",
    "LINEARVISCOELASTICORTHOTROPICPOWERLAW",
    "LINEARELASTIC",
    "LINEARVISCOELASTICPOWERLAW",
    "CDP",
    "B4",
    "ADVONMISES",
    "ADLINEARELASTIC",
]
finitestrain_materials = [
    "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY",
    "FINITESTRAINJ2PLASTICITY",
    "FINITESTRAINISOTROPICBIOTVISCOELASTICITY",
    "FINITESTRAINORTHOTROPICBIOTVISCOELASTICITY_SUBSTEPPED",
    "COMPRESSIBLENEOHOOKE",
    "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY_SUBSTEPPED",
    "FINITESTRAINJ2PLASTICITY_SUBSTEPPED",
    "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY",
    "ADCOMPRESSIBLENEOHOOKE",
]

material_dirs = os.listdir("modules/materials")


def find_props(mat_dir):
    test_path = os.path.join("modules/materials", mat_dir, "test", "test.cpp")
    if not os.path.exists(test_path):
        return None
    with open(test_path, "r") as f:
        content = f.read()

    # Looking for materialProperties = { ... } or double props[] = { ... }
    m = re.search(r"materialProperties\s*(?:\[\])?\s*=\s*\{([^\}]+)\}", content)
    if m:
        return m.group(1).strip()

    m = re.search(r"double\s+props(?:\[\])?\s*=\s*\{([^\}]+)\}", content)
    if m:
        return m.group(1).strip()

    m = re.search(r"std::vector<\s*double\s*>\s+\w+\s*=\s*\{([^\}]+)\}", content)
    if m:
        return m.group(1).strip()

    return None


os.makedirs("examples/python", exist_ok=True)

for mat_name in hypoelastic_materials:
    # Find matching directory
    mat_dir = None
    for d in material_dirs:
        if d.upper() == mat_name or d.upper() + "MODEL" == mat_name or d.upper() == mat_name.replace("_SUBSTEPPED", ""):
            mat_dir = d
            break

    props = "1.0, 0.3"  # default fallback
    if mat_dir:
        p = find_props(mat_dir)
        if p:
            props = p

    code = f"""import sys
import os
import numpy as np

# Add the build lib directory to PYTHONPATH to find marmot extension
sys.path.insert(0, os.path.abspath('../../build/lib'))

import marmot

print("Running example for HypoElastic material: {mat_name}")

# Material properties extracted from C++ tests
properties = np.array([{props}], dtype=np.float64)

# Setup solver
options = marmot.solvers.HypoElasticSolver.SolverOptions()
solver = marmot.solvers.HypoElasticSolver("{mat_name}", properties, options)

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
    print(f"Final Stress 11: {{history[-1].stress[0]}}")
"""
    with open(f"examples/python/example_{mat_name.lower()}.py", "w") as f:
        f.write(code)


for mat_name in finitestrain_materials:
    # Find matching directory
    mat_dir = None
    for d in material_dirs:
        if d.upper() == mat_name or d.upper() + "MODEL" == mat_name or d.upper() == mat_name.replace("_SUBSTEPPED", ""):
            mat_dir = d
            break

    props = "1.0, 0.3"  # default fallback
    if mat_dir:
        p = find_props(mat_dir)
        if p:
            props = p

    code = f"""import sys
import os
import numpy as np

# Add the build lib directory to PYTHONPATH to find marmot extension
sys.path.insert(0, os.path.abspath('../../build/lib'))

import marmot

print("Running example for FiniteStrain material: {mat_name}")

# Material properties extracted from C++ tests
properties = np.array([{props}], dtype=np.float64)

# Setup solver
options = marmot.solvers.FiniteStrainSolver.SolverOptions()
solver = marmot.solvers.FiniteStrainSolver("{mat_name}", properties, options)

# Setup a loading step
step = marmot.solvers.FiniteStrainSolver.Step()
step.gradUIncrementTarget = np.array([
    [0.01, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0]
], dtype=np.float64)

step.isGradUComponentControlled = np.array([
    [True,  True,  True],
    [True,  True, True],
    [True,  True,  True]
])
step.isStressComponentControlled = np.logical_not(step.isGradUComponentControlled)
step.stressIncrementTarget = np.zeros((3,3), dtype=np.float64)

solver.addStep(step)
solver.solve()

history = solver.getHistory()
if len(history) > 0:
    print(f"Final Stress XX: {{history[-1].stress[0,0]}}")
"""
    with open(f"examples/python/example_{mat_name.lower()}.py", "w") as f:
        f.write(code)

print("Generated all python examples in examples/python/")
