"""
3,6,9j symbols calculator 
Filename: wigner_symbols.py
Author: Neville Chen
Date: 01 Jul 2021

Script to calculate 6j symbols.
Explicit formula taken from Dan Steck.

"""

import numpy as np 
from math import factorial
from . import cg_coeff 

def sum_range_6j(j,k1,k2,k3,m1,m2,m3):
	"""
	Define range for summation 
	"""
	start = int(max([j,k1,k2,k3]))
	end = int(min([m1,m2,m3]))+1 #range function excludes endpoint
	return range(start,end)

def delta(a,b,c):
	if a+b-c<0 or b+c-a<0 or c+a-b<0:
		return 0
	else:
		num = factorial(a+b-c)*factorial(b+c-a)*factorial(c+a-b)
		denom = factorial(a+b+c+1)
		return np.sqrt(num/denom)

def wigner_3j(j1,m1,j2,m2,j3,m3):
	"""
	Evaluate 3j symbols according to
	j1 j2 j3
	m1 m2 m3
	"""
	return (-1)**(j1-j2-m3)/np.sqrt(2*j3+1)*cg_coeff.cg(j1,m1,j2,m2,j3,-m3)


def wigner_6j(j1,j2,j3,l1,l2,l3):
	"""
	Wigner 6j evaluated according to
	j1 j2 j3
	l1 l2 l3

	The mapping to actual momenta is 
	j1 j2 j12
	j3 j  j23
	"""
	j=j1+j2+j3
	k1=j1+l2+l3
	k2=l1+j2+l3
	k3=l1+l2+j3
	m1=j1+j2+l1+l2
	m2=j2+j3+l2+l3
	m3=j3+j1+l3+l1

	to_mul = delta(j1,j2,j3)*delta(j1,l2,l3)*delta(l1,j2,l3)*delta(l1,l2,j3)

	srange = sum_range_6j(j,k1,k2,k3,m1,m2,m3)

	val=0
	for n in srange:
		val+=(-1)**n*factorial(n+1)/(factorial(n-j)*factorial(n-k1)*factorial(n-k2)*factorial(n-k3)*factorial(m1-n)*factorial(m2-n)*factorial(m3-n))

	return to_mul*val

def recoup3(j1, j2, j3, j12, j23, j):
	"""
	gives \braket{j12j3j}{j1j23j}
	"""
	return (-1)**(j1+j2+j3+j)*np.sqrt((2*j12+1)*(2*j23+1))*wigner_6j(j1,j2,j12,j3,j,j23)




