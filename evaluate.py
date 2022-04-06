nsites = 500

'''
values = {}

tot = 0

for k in range(2 << (nsites - 1)):
	config = [1 if v == '1' else -1 for v in bin(k)[2:].zfill(nsites)]
	f_sum = sum(config[i - 1] * config[i] for i in range(len(config)))
	s_sum = sum(config)
	values[(s_sum, f_sum)] = []
for k in range(2 << (nsites - 1)):
	#s = SpinConfig(bin(k)[2:].zfill(nsites))
	#tot += (np.exp(-1 * H.energy(s)))
	config = [1 if v == '1' else -1 for v in bin(k)[2:].zfill(nsites)]
	f_sum = sum(config[i - 1] * config[i] for i in range(len(config)))
	s_sum = sum(config)
	values[(s_sum, f_sum)].append(bin(k)[2:].zfill(nsites))

for k in sorted(values):
	print(k, len(values[k]))

#print(tot)
'''

from math import comb
import numpy as np


J = -1
mu = 1.1
tot = 0

n = 0
for sec_sum in range(nsites, -1, -2):
	first_sum = nsites - 4
	#print("N:", n)
	if n == 0:
		energy_neg = np.exp(-1 * (-J * nsites + mu * sec_sum))
		#print(sec_sum, nsites, 1)
		energy_pos = np.exp(-1 * (-J * nsites - mu * sec_sum))
		#print(-1 * sec_sum, nsites, 1)
		tot += (energy_pos + energy_neg)
	elif n == 1:
		energy_pos = np.exp(-1 * (-J * first_sum + mu * sec_sum))
		#print(sec_sum, first_sum, nsites)
		energy_neg = np.exp(-1 * (-J * first_sum - mu * sec_sum))
		#print(-1 * sec_sum, first_sum, nsites)
		tot += (nsites * energy_pos + nsites * energy_neg)
	else:
		num_of_zrs = nsites - n
		for k in range(1, n + 1):
			freq = int((comb(num_of_zrs - 1, k - 1) * comb(n - 1, k - 1) * nsites)/k)
			energy_pos = np.exp(-1 * (-J * first_sum + mu * sec_sum))
			tot += freq * energy_pos
			#print(sec_sum, first_sum, freq)
			if n != nsites//2 or nsites % 2 != 0:
				energy_neg = np.exp(-1 * (-J * first_sum - mu * sec_sum))
				tot += freq * energy_neg
				#print(-1 * sec_sum, first_sum, freq)

			first_sum -= 4
	n += 1

print(tot)

'''
from montecarlo import *

H = IsingHamiltonian(k=1, J=-1, mu=1.1)
H.compute_average_energy(SpinConfig("+" * nsites), 1)
'''