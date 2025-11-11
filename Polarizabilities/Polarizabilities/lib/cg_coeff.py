"""
Clebsch-Gordan coefficients and 3j symbols calculator
Filename: cg_coeff.py
Author: Neville Chen
Date: 4 June 2021

Script to calculate Clebsch-Gordan coefficients using the explicit formula taken from Dan Steck's notes (Eqn 7.51).

Inputs:
j1 - Scalar for first total angular momentum
m1 - Scalar for projection of j1
j2 - Scalar for second total angular momentum
m2 - Scalar for projection of j2
j3 - Scalar for total coupled angular momentum
m3 - Scalar for projection of j3

Outputs:
cg - Scalar for the Clebsch Gordan coefficient. It includes the sign and square root.

"""


import numpy as np 
from math import factorial 


##################### Auxiliary functions ###########################

def sum_range_cg(j1,m1,j2,m2,j3,m3):
	"""
	Define range for summation 
	"""
	start = int(max([0,j2-j3-m1,j1+m2-j3]))
	#range function excludes endpoint
	end = int(min([j1-m1,j2+m2,j1+j2-j3]))+1 
	return range(start,end)

def kronecker(m1,m2,m3):
	"""
	Delta function to ensure that m1+m2=m3
	"""
	return int((m1+m2) != m3)

def triangle(j1,j2,j3):
	"""
	Check that j1,j2,j3 fulfil the triangle relation
	"""
	list_j3 = np.arange(np.abs(j1-j2),j1+j2+1,1)
	return int(j3 not in list_j3)

def delta(a,b,c):
	if a+b-c<0 or b+c-a<0 or c+a-b<0:
		return 0
	num = factorial(a+b-c)*factorial(b+c-a)*factorial(c+a-b)
	denom = factorial(a+b+c+1)
	return np.sqrt(num/denom)

def check_valid(j1,m1,j2,m2,j3,m3):
	"""
	Function to check that choices of j-s and m-s are valid

	1. j-s are all >=0
	2. j-s are half or full integers
	3. m-s are bounded between -j to j
	4. m-s are half or full integers
	"""
	flag_pos_j = j1<0 or j2<0 or j3<0
	flag_int_j = (4*j1)%2 != 0 or (4*j2)%2 != 0 or (4*j3)%2 != 0
	flag_mag_m = np.abs(m1)>j1 or np.abs(m2)>j2 or np.abs(m3)>j3
	flag_int_m = (4*m1)%2 != 0 or (4*m2)%2 != 0 or (4*m3)%2 != 0
	flag_kronecker = kronecker(m1,m2,m3)
	flag_tri = triangle(j1,j2,j3)

	flag = flag_pos_j + flag_mag_m + flag_int_j +flag_int_m +flag_kronecker+flag_tri
	

	# if flag_pos_j == True:
	# 	print("All j-s must be positive")
	# if flag_int_j == True:
	# 	print("j-s must be half or full integer")
	# if flag_mag_m == True:
	# 	print("All m-s must be smaller than their respective j-s")
	# if flag_int_m == True:
	# 	print("All m-s must be full or half integer")
	# if flag_kronecker == True:
	# 	print("m1+m2 is not equal to m3")
	# if flag_tri == True:
	# 	print("j3 is not a valid result of addition")

	if flag != 0:
		return 0

	else:
		return 1

#########################################################################


def cg(j1,m1,j2,m2,j3,m3):
	"""
	Main function to calculate Clebsch Gordan coefficients
	"""
	flag_valid = check_valid(j1,m1,j2,m2,j3,m3)
	if flag_valid == 0:
		# print("Please check the above error(s)")
		return 0


	#convert to float so that np.sqrt can handle large numbers
	large_factorial=float(factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j3+m3)*factorial(j3-m3)*(2*j3+1)) 
	to_mul= delta(j1,j2,j3)*np.sqrt(large_factorial)

	#determine range of sum
	srange = sum_range_cg(j1,m1,j2,m2,j3,m3)
	
	#evaluate sum
	val = 0
	for k in srange:
		val+=(-1)**k/(factorial(j1-m1-k)*factorial(j3-j2+m1+k)*factorial(j2+m2-k)*factorial(j3-j1-m2+k)*factorial(j1+j2-j3-k)*factorial(k))

	return to_mul*val



