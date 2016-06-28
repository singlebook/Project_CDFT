#!/usr/bin/python
#Read the parameters of the system from a file

import sys 
from phys_para import *

class SystemDef:
	'''
	approximation for attractive excess helmholtz free energy
    MSA:  Mean Spherical Approximation
    FMSA: First order Mean Spherical Approximation
    MFA: Mean field approximation


    approximation for estimating effective HS size
    BH :   Baker-Handerson method
    WCA:   Weitheim-Chandler-Anderson method
    OCB:   Oxtoby



    reasonable combination:
    BH  + MFA
    BH  + FMSA
    '''
	temperature = 0
	beta = 0 # 1 / kT
	rcut = 0  # cutoff for LJ attraction
	size = [] # The size of the full system, size[0] is X dirextion, size[1] is Y dirextion, size[2] is Z dirextion.
	number_component = 0
	bulk_density = []
	LJ_parameters = [] # [epsion_0, sigma_0, epsion_1, sigma_1,.....]
	pts = [] # Discretize the space into a given number of grid points, pts[0] is X dirextion, pts[1] is Y dirextion, pts[2] is Z dirextion
	delta = [] #The distance between two grid nodes
	atom_type = ""
	HS_size_approximation = ""
	approx_excess_free_energy = ""
	input_external_field = ""
	
	def __init__(self, name):
		self.name = name   
		input = open(self.name)
		for line in input:
			paras = line.split()
			if paras[0] == "Temperature":
				self.temperature = float(paras[1])
				self.beta = 1.0 / (const.kb * const.Na * self.temperature)
			elif paras[0] == "Size":
				self.size = [float(paras[1]), float(paras[2]), float(paras[3])]
			elif paras[0] == "Lattice":
				self.pts = [float(paras[1]), float(paras[2]), float(paras[3])]
				self.delta = [self.size[0]/(self.pts[0]-1), self.size[1]/(self.pts[1]-1), self.size[2]/(self.pts[2]-1)]
			elif paras[0] == "AtomType":
				self.atom_type = paras[1]
	   		elif paras[0] == "SizeApproximation":
				self.HS_size_approximation = paras[1]
			elif paras[0] == "CutOff":
				self.rcut = float(paras[1])
			elif paras[0] == "Attractive_excess_helmholtz_free_energy":
			    self.approx_excess_free_energy = paras[1]
			elif paras[0] == "Component":
				self.number_component = int(paras[1])
			elif paras[0] == "Bulk_Density":
				if self.number_component <= 0:
					print "Error:  The number of components must be positive and defined before the defination of bulk densities. \n"
					sys.exit() 
				elif self.number_component > 0:
					self.bulk_density = paras[1:self.number_component+1]
					for i in range(len(self.bulk_density)):
						self.bulk_density[i] = float(self.bulk_density[i])
			elif paras[0] == "LJ_Para":
				if self.number_component <= 0:
					print "Error:  The number of components must be positive and defined before the defination of LJ parameters. \n"
					sys.exit() 
				elif self.number_component > 0:
					self.LJ_parameters = paras[1:2*self.number_component+1]
					for i in range(len(self.LJ_parameters)):
						self.LJ_parameters[i] = float(self.LJ_parameters[i])
			elif paras[0] == "External_Field":
				self.input_external_field = paras[1]
				
				
			    
