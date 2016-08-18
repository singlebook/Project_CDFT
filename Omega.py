#!/usr/bin/python
#This module contains the Fourier transforms of the weight functions.
#Author: Wei Chen
#Date: Aug 18, 2016
#Place : Beijing 

import math

def Omega0(K, R):
	return math.sin(K*R) / (K*R)
	
def Omega1(K, R):
	return math.sin(K*R) / K

def Omega2(K, R):
	return 4.0 * math.pi * R * math.sin(K*R)/K
	
def Omega3(K, R):
	return (4.0 * math.pi / (K*K*K)) * (math.sin(K*R) - K*R*math.cos(K*R))
	
