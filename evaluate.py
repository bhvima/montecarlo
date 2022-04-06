nsites = 10

values = {}

for k in range(2 << (nsites - 1)):
	config = [1 if v == '1' else -1 for v in bin(k)[2:].zfill(nsites)]
	f_sum = sum(config[i - 1] * config[i] for i in range(len(config)))
	s_sum = sum(config)
	values[(s_sum, f_sum)] = []
for k in range(2 << (nsites - 1)):
	config = [1 if v == '1' else -1 for v in bin(k)[2:].zfill(nsites)]
	f_sum = sum(config[i - 1] * config[i] for i in range(len(config)))
	s_sum = sum(config)
	values[(s_sum, f_sum)].append(bin(k)[2:].zfill(nsites))

for k in sorted(values):
	print(k, len(values[k]))

#print(len(values))

from math import comb

n = 0
for sec_sum in range(nsites, -1, -2):
	first_sum = nsites - 4
	print("N:", n)
	if n == 0:
		print(1, (sec_sum, nsites))
	elif n == 1:
		print(nsites, (sec_sum, first_sum))
	else:
		num_of_zrs = nsites - n
		for k in range(1, n + 1):
			freq = int((comb(num_of_zrs - 1, k - 1) * comb(n - 1, k - 1) * nsites)/k)
			print(freq, (sec_sum, first_sum))
			first_sum -= 4
	n += 1