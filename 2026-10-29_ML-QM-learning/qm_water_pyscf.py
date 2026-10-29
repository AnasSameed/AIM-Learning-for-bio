"""
DFT practice: optimize water, then confirm it is a minimum (no imaginary frequencies).
B3LYP / def2-SVP using PySCF.   pip install pyscf geometric
Expected: 3 real vibrations near 1600, 3700, 3800 cm-1 (values depend on the basis set).
NOTE: written as a learning script; run it in Colab and check the output yourself.
"""
import numpy as np
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize
from pyscf.hessian import thermo

mol = gto.M(atom="O 0 0 0; H 0 0.76 0.59; H 0 -0.76 0.59", basis="def2-svp", verbose=0)
mf = dft.RKS(mol); mf.xc = "b3lyp"
mol_eq = optimize(mf)                                   # 1. geometry optimization

mf2 = dft.RKS(mol_eq); mf2.xc = "b3lyp"
e = mf2.kernel()                                        # 2. single-point energy
print(f"Energy (Hartree): {e:.6f}")

hess = mf2.Hessian().kernel()                           # 3. frequency calculation
freqs = thermo.harmonic_analysis(mol_eq, hess)["freq_wavenumber"]
print("Frequencies (cm-1):", np.round(freqs, 1))
n_imag = int(np.sum(np.imag(freqs) != 0))
print("Imaginary frequencies:", n_imag, "-> minimum" if n_imag == 0 else "-> saddle point")
# Transition state check would instead expect exactly ONE imaginary frequency, then an IRC.
