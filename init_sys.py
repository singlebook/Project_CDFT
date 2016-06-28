#!/usr/bin/python
#Initiate the external field and the density profile
#Author: Wei Chen
#Date: June 28, 2016
#Place : Beijing

import sys
import system_def
import math
from phys_para import *

class Init():
	field_external = []
	density = []
	def __init__(self, filename, parameters):
		self.name = filename
		self.parameters = parameters
		input = open(self.name)
		for line in input:
			paras = line.split()
			self.field_external.append(float(paras[3]))
			self.density.append(math.exp(-1.0 * self.parameters * float(paras[3])))
		input.close()
			
