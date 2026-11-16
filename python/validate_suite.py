from pathlib import Path
import numpy as np

from solver_suite import (
    run_cpo_dbo, run_cpo, run_dbo, run_ga, run_ipso,
    run_ga_pso, run_aco_ga, fitness,
)


ROOT = Path(__file__).resolve().parents[1] / "data"
dose = np.loadtxt(ROOT / "instance_lib20" / "doseMatrix_inst01.csv", delimiter=",")

tests = [
    ("CPO-DBO", lambda: run_cpo_dbo(dose, 1234, m=40, budget=8_000)),
    ("CPO", lambda: run_cpo(dose, 1234, m=40, budget=8_000)),
    ("DBO", lambda: run_dbo(dose, 1234, m=40, budget=8_000)),
    ("GA", lambda: run_ga(dose, 1234, m=40, budget=8_000)),
    ("IPSO", lambda: run_ipso(dose, 1234, m=40, budget=8_000)),
    ("GA-PSO", lambda: run_ga_pso(dose, 1234, m=40, budget=8_000)),
    ("ACO-GA", lambda: run_aco_ga(dose, 1234, m=40, budget=8_000)),
]

hashes = {}
for name, call in tests:
    result = call()
    value, path, curve, first_hit, evaluations = result[:5]
    initial = result[-1]
    assert np.array_equal(np.sort(path), np.arange(1, 51)), name
    assert abs(fitness(path, dose) - value) < 1e-9, name
    assert np.all(np.diff(curve) <= 1e-10), name
    assert 1 <= first_hit <= 10_000, name
    hashes[name] = hash(initial.tobytes())
    print(name, round(float(value), 6), evaluations, first_hit)

common = [hashes[k] for k in ["CPO-DBO", "CPO", "DBO", "GA", "GA-PSO", "ACO-GA"]]
assert len(set(common)) == 1, hashes
print("VALIDATION PASSED")
