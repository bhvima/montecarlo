import timeit
from montecarlo import *

H = IsingHamiltonian(-1, -2, 1.1)

for i in range(1, 30, 1):
	s = SpinConfig("+" * i)
	print(i, timeit.timeit("H.compute_average_energy(s, 2)", globals=globals(), number=1))


