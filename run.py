#!/usr/bin/python
import bkevol as bk

S = bk.Settings()

#Evolution settings
S.alphaS = 0.118
S.freezingScale = 1.0
S.eps  = 1e-7
S.mu2  = 1e-2


S.NkT2 = 17
S.NkT2int = 17

#S.NkT2 = 257
#S.NkT2int = 257


S.kT2Min = 1e-2
S.kT2Max = 1e6


S.Nrap = 64
#S.Nrap = 512
S.xMin = 1e-6
S.xMax = 1



#S.kernelType = "BFKLplain:sub" #)
#S.kernelType = "BFKL:eps" #+sub)
S.kernelType = "BFKL_res:Eps" #+sub)
#S.kernelType = "BFKL_res_kc_simp:eps" #+sub)
#S.kernelType = "BFKL_res_kc_v_r_simp:eps" #+sub)
#S.kernelType = "BFKL_res_kc_full:eps" #+sub)
#S.kernelType = "BFKL_res_kc_v_r_full:eps" #+sub)
#S.kernelType = "BFKL_res_DGLAP:eps"
#S.kernelType = "BFKL_res_kc_full_DGLAP:eps"

S.funStr = "function"
S.pars = [(1,2, 3),
          (2,4., 5)]


#S.printInfo()
#import sys
#sys.exit(0)

kernels = [
 "BFKLplain:sub",
 "BFKL:Eps", "BFKL:sub",
 "BFKL_res:Eps", "BFKL_res:sub",
 "BFKL_res_kc_simp:Eps", "BFKL_res_kc_simp:sub",
 "BFKL_res_kc_v_r_simp:Eps", "BFKL_res_kc_v_r_simp:sub",
 "BFKL_res_kc_full:Eps", "BFKL_res_kc_full:sub",
 "BFKL_res_kc_v_r_full:Eps", "BFKL_res_kc_v_r_full:sub",
 "BFKL_res_DGLAP:Eps",
 "BFKL_res_kc_full_DGLAP:Eps",
 "BFKL_res_kc_full_DGLAP_simp_kc:Eps",
 "BFKL_res_kc_full_DGLAP_full_kc:Eps",
]


#S.kernelType = "BFKLplain:Sub"
S.kernelType = "BFKL_res_DGLAP:Eps"
#S.kernelType = "BFKL_res_DGLAP:ZEps"

#S.kernelType = "BFKL_res_kc_full_DGLAP:Eps"
#S.kernelType = "BFKL_res_kc_full_DGLAP_simp_kc:Eps"
#S.kernelType = "BFKL_res_kc_full_DGLAP_simp_kc:ZEps"

#S.kernelType = "BFKL_res_kc_full_DGLAP_full_kc:Eps"

print S.alphaS

#for k in kernels:
#S.kernelType = k
sol = bk.Solver(S)
sol.CalcEvolKernel()

#sol.EvolveAll()
#sol.PrintBaseGrid();

sol.SaveEvolKernels("testNow.h5")
#sol.LoadEvolKernels("testNow.h5")
#sol.EvolveAll()
sol.PrintBaseGrid();
