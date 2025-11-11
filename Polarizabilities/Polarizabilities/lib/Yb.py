import json
import numpy as np
from pathlib import Path

# ensures that the path for YbAtomicData works for all users
file_dir = Path( __file__ ).parent
data_path = file_dir.joinpath("YbAtomicData.json")

class Yb:

	#  constants for atomic structure
	Z = 70
	# following convention of Hz = gS * S + gL * L - gI * I
	gS = 2.0023193043737 # g-factor of electronic spin
	gL = 1.0 # g-factor of orbital angular momentum 
	gI = 0.0 #g-factor of nuclear spin

	# isotope data
	isotope_abundance = {
	"171": 0.143,
	"174": 0.318
	}

	isotope_mass = {
	"171": 170.936323,
	"174": 173.938859
	}

	def __init__(self):
		self.data = json.load(open(data_path))
		self.ground_state = 6 

	# creating the key for a certain orbital
	def encode(self, n, l, s, j, type_str='o'): 
		return type_str + str(n) + str(l) + str(s) + str(j)

	# decode key to find quantum numbers
	def decode(self, key_str):
		n = int(key_str[1])
		l = int(key_str[2])
		s = int(key_str[3])
		j = int(key_str[4])
		return (n , l, s, j, key_str[0])

	def getStateFromQnumbers(self, n, l, s, j):
		L_symbol = ["S", "P", "D", "F"]
		str_L = L_symbol[l]
		str_S = str(int(2*s +1))
		str_J = str(j)
		str_n = "n=" + str(n) +", "
		return str_n + str_S + str_L + str_J

	def getStateFromKey(self, key):
		return self.getStateFromQnumbers(*self.decode(key))
	
	#returns energy of a certain orbital 
	def getEnergy(self, n, l, s, j, type_str='o'):
		if n < self.ground_state:
			raise ValueError("n must be >= 6")
		if l >= n:
			raise ValueError("l must not be greater than n")

		try:
			return self.data[self.encode(n,l,s,j, type_str)]["Energy"]
		except KeyError:
			raise KeyError("Invalid state")
			return None

	#returns the RDME of the 1st orbital relative to the second 
	def getRDME(self, n1, l1, s1, j1, n2, l2, s2, j2, type_str2):
		if n1 < self.ground_state:
			raise ValueError("n1 must be >= 6")
		if l1 >= n1:
			raise ValueError("l must not be greater than n")
		if l2 >= n2:
			raise ValueError("l2 must not be greater than n2")

		try:
			return self.data[self.encode(n1,l1,s1,j1)]["RDME"][self.encode(n2,l2,s2,j2,type_str2)]
		except KeyError:
			print("RDME for state not in database")
			return None

	def getRDMEKeys(self, n, l, s, j, type_str='o'):
		if n < self.ground_state:
			raise ValueError("n must be >= 6")
		return [*self.data[self.encode(n,s,l,j,type_str)]["RDME"]]

	#Adds a new orbital with an associated energy value by overwriting the old database
	def addEnergy(self, n, l, s, j, energy):
		if n < 6:
			raise ValueError("n must be >= 6")
		if l >= n:
			raise ValueError("l must not be greater than n")
		temp = {
			"Energy": energy,
		}
		self.data.update( {self.encode(n,l,s,j) : temp} )

		j = json.dumps(self.data, indent=2) # indent is for readability
		with open(data_path, 'w') as f:
			f.write(j)
			f.close()

	#Adds another RDME to the 1st orbital by overwriting the old database
	def addRDME(self, n1, l1, s1, j1, n2, l2, s2, j2, rdme):
		if n1 < 6:
			raise ValueError("n1 must be >= 6")
		if l1 >= n1:
			raise ValueError("l must not be greater than n")
		if l2 >= n2:
			raise ValueError("l2 must not be greater than n2")

		# add RDME to first state
		if "RDME" in self.data[self.encode(n1,s1,l1,j1)]:
			self.data[self.encode(n1,l1,s1,j1)]["RDME"].update({self.encode(n2,l2,s2,j2) : rdme})
		else:
			temp ={ "RDME" : {self.encode(n2,l2,s2,j2) : rdme}}
			self.data[self.encode(n1,l1,s1,j1)].update(temp)

		# add RDME to second state
		if "RDME" in self.data[self.encode(n2,l2,s2,j2)]:
			self.data[self.encode(n2,l2,s2,j2)]["RDME"].update({self.encode(n1,l1,s1,j1) : rdme})
		else:
			temp ={ "RDME" : {self.encode(n1,l1,s1,j1) : rdme}}
			self.data[self.encode(n2,l2,s2,j2)].update(temp)


		j = json.dumps(self.data, indent=2) # indent is for readability
		with open(data_path, 'w') as f:
			f.write(j)
			f.close()


# atom = Yb()
# print(atom.getStateFromQnumbers(6,2,1,1))
# atom.addRDME(6,1,1,2,7,0,1,1,5.05)