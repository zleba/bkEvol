#!/usr/bin/python
import bkevol as bk

S = bk.Settings()

#S.asMZ = 0.2
#S.eps  = 1e-7

S.N = 33
S.Nint = 33
S.Nrap = 64

S.inputDir = "BFKLplain:sub"

print S.asMZ

sol = bk.Solver(S)
sol.InitMat()
sol.EvolveNew()
sol.PrintBaseGrid();
