"""Local algebra/source tests; these do NOT replace a Cactus/MPI build."""
import ast
import re
import subprocess
import tempfile
import unittest
from pathlib import Path
import os

import numpy as np
from em_constraint_norms import normalized_rows

THORN = Path(__file__).resolve().parents[1]
SRC = THORN / "src"


def uncomment(text):
    return re.sub(r"\s+", "", "\n".join(line.split("!", 1)[0] for line in text.splitlines()))


def gauss_python():
    lines, indent = [], 0
    arrays = ("cf1", "cf2", "hu", "d1_hh", "d1_lE", "d1_lB", "lE", "lB", "d1_ch")
    for raw in (SRC / "EinsteinMaxwell_Gauss.inc").read_text().splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if line == "end do":
            indent -= 1
            continue
        match = re.fullmatch(r"do (\w+) = (\w+), (\w+)", line)
        if match:
            var, start, end = match.groups()
            lines.append("    " * indent + f"for {var} in range({start}, {end}+1):")
            indent += 1
            continue
        line = re.sub(r"(\d+(?:\.\d*)?)d0", r"\1", line)
        for name in arrays:
            line = re.sub(rf"\b{name}\(([^()]*)\)", rf"{name}[\1]", line)
        line = re.sub(r"^(cf[12]) =", r"\1[:] =", line)
        lines.append("    " * indent + line)
    return compile("\n".join(lines), "Gauss.inc algebra", "exec")


class DiagnosticsTests(unittest.TestCase):
    def test_rhs_numerics_unchanged(self):
        relative = "einstein-maxwell/EinsteinMaxwellEvolve/src/EinsteinMaxwell_calc_rhs.F90"
        reference = os.environ.get(
            "DYON_RHS_REFERENCE",
            "0dabbe9c9e3b5a2524fd259dd2a2f7d7578058d2",
        )
        original = subprocess.check_output(
            ["git", "show", f"{reference}:{relative}"],
            cwd=THORN,
            text=True,
        )
        current = (SRC / "EinsteinMaxwell_calc_rhs.F90").read_text()
        for name in ("EinsteinMaxwell_FD.h", "EinsteinMaxwell_Gauss.inc"):
            current = current.replace(f'#include "{name}"', (SRC / name).read_text())
        self.assertEqual(uncomment(original), uncomment(current))

    def test_diagnostic_derivatives_match_rhs(self):
        pattern = r"^\s*(d1_(?:ch|hh|lE|lB)\([^\n]*\))\s*=\s*([^\n]+)"
        rhs = re.findall(pattern, (SRC / "EinsteinMaxwell_calc_rhs.F90").read_text(), re.M)
        diag = re.findall(pattern, (SRC / "EinsteinMaxwell_calc_constraints.F90").read_text(), re.M)
        self.assertEqual([(uncomment(a), uncomment(b)) for a,b in rhs],
                         [(uncomment(a), uncomment(b)) for a,b in diag])

    def test_gauss_against_metric_determinant_identity(self):
        rng = np.random.default_rng(17)
        code = gauss_python()
        for _ in range(100):
            q = rng.normal(size=(3, 3))
            h = q @ q.T + np.eye(3)
            dh = rng.normal(size=(3, 3, 3))
            dh = (dh + dh.swapaxes(0, 1)) / 2
            env = {name: np.zeros(shape) for name, shape in (
                ("cf1",(4,4,4)),("cf2",(4,4,4)),("hu",(4,4)),
                ("d1_hh",(4,4,4)),("d1_lE",(4,4)),("d1_lB",(4,4)),
                ("lE",(4,)),("lB",(4,)),("d1_ch",(4,)))}
            env["hu"][1:,1:] = np.linalg.inv(h)
            env["d1_hh"][1:,1:,1:] = dh
            for name in ("lE", "lB", "d1_ch"):
                env[name][1:] = rng.normal(size=3)
            for name in ("d1_lE", "d1_lB"):
                env[name][1:,1:] = rng.normal(size=(3,3))
            env.update(ch=0.7, conf_fac_exponent=1.0)
            connection = np.array([0.5*np.trace(np.linalg.inv(h) @ dh[:,:,k]) for k in range(3)])
            connection -= 1.5 * env["d1_ch"][1:] / env["ch"]
            expected = [np.trace(env[d][1:,1:]) + connection @ env[v][1:]
                        for d,v in (("d1_lE","lE"),("d1_lB","lB"))]
            exec(code, env)
            np.testing.assert_allclose([env["divE"],env["divB"]], expected, rtol=1e-12, atol=1e-12)

    def test_grouped_scalar_columns_and_mask_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for reduction in ("norm2", "norm_inf"):
                for kind, values in (("e", "2 99"), ("b", "3 99")):
                    (root / f"group_{kind}.{reduction}.asc").write_text(
                        f"# data columns: 3:EM_C{kind.upper()}_near 4:EM_C{kind.upper()}_far\n0 0 {values}\n")
            (root / "group_mask.average.asc").write_text(
                "# data columns: 3:EM_mask_near 4:EM_mask_far\n0 0 0.25 0.75\n")
            rows = list(normalized_rows(root, "near"))
            self.assertEqual(rows[0][2:6], (4.0,6.0,2.0,3.0))

    def test_generator_syntax(self):
        generator = os.environ.get("DYON_GENERATOR")
        if not generator:
            self.fail("Set DYON_GENERATOR to the generator path.")
        ast.parse(Path(generator).read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
