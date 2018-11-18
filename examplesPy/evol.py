#!/usr/bin/env python
import bkevol as bk

# Load Setting object with configuration of the evolution

S = bk.Settings()

# Value of the $\alpha_S(M_Z)$ and the $\alpha_S$ freezing scale

S.alphaS = 0.118
S.freezingScale = 1.0

# The regularization variable eps and the scale $\mu$

S.eps  = 1e-7
S.mu2  = 1e-2

# Number of the grid points in $k_T^2$ variable and the range in $k_T^2$

S.NkT2 = S.NkT2int = 17
S.kT2Min = 1e-2
S.kT2Max = 1e6

# Number of the grid points in $x$ variable and the xMin and xMax values

S.Nrap = 64
#S.Nrap = 512
S.xMin = 1e-6
S.xMax = 1

# Function to be evolved (TODO)

S.funStr = "function"
S.pars = [(1,2, 3),
          (2,4., 5)]

# Overview of the implemented kernels

kernels = [
 "BFKLplain:Sub",
 "BFKL:Eps", "BFKL:Sub",
 "BFKL_res:Eps", "BFKL_res:Sub",
 "BFKL_res_kc_simp:Eps", "BFKL_res_kc_simp:Sub",
 "BFKL_res_kc_v_r_simp:Eps", "BFKL_res_kc_v_r_simp:Sub",
 "BFKL_res_kc_full:Eps", "BFKL_res_kc_full:Sub",
 "BFKL_res_kc_v_r_full:Eps", "BFKL_res_kc_v_r_full:Sub",
 "BFKL_res_DGLAP:Eps",
 "BFKL_res_kc_full_DGLAP:Eps",
 "BFKL_res_kc_full_DGLAP_simp_kc:Eps",
 "BFKL_res_kc_full_DGLAP_full_kc:Eps",
]

# Evolution kernel to be used

S.kernelType = "BFKL:Sub"

# Calculate the evolution kernel

sol = bk.Solver(S)
sol.CalcEvolKernel()

# Posibility to save the evolution kernel

sol.SaveEvolKernels("evolKernel.h5")
#sol.LoadEvolKernels("evolKernel.h5")

# Evolve with the function in the Settings

sol.EvolveAll()

# Get the result of the evolution and the grid points

#sol.PrintBaseGrid()
res = sol.getResult()
k2Grid = sol.getK2grid()

# Plot results

#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
for i,r in enumerate(res[::6]):
    plt.loglog(k2Grid, r, label=i)
plt.ylim(1e-7, 1e6)
plt.xlabel('$k^2 [GeV^2]$')
plt.ylabel('$\phi$')
plt.show()
