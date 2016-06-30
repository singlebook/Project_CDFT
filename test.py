#!/usr/bin/python
#testing

import sys
import system_def
import init_sys
import free_energy_func
from scipy.optimize import minimize
import math

para = system_def.SystemDef(sys.argv[1])

init = init_sys.Init(para.input_external_field, para.beta)

field_ext =init.field_external

init.density

def F_id_ext(x):
	fid = 0
	fext = 0
	for i in range(0, para.pts[2]+1):
			for j in range(0, para.pts[1]+1):
				for k in range(0, para.pts[0]+1):
					index = k + para.pts[0] * (j + i * para.pts[1])
					fid += para.delta[0] * para.delta[1] * para.delta[2] * x[index] * (math.log(x[index]) - 1.0)   
					fext += para.delta[0] * para.delta[1] * para.delta[2] * x[index] * field_ext[index] 
	return (fid+fext)
					
res = minimize(F_id_ext, init.density, method='L-BFGS-B', jac = False, tol=1e-4)					 
#f = free_energy_func.F_id_ext(para.beta, para.pts, density, field_ext, para.delta, para.thermal_wave)


print "The simulation is finally over!!"
