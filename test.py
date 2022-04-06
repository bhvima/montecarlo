from montecarlo import *

N = 500

H = IsingHamiltonian(k=1, J=-1, mu=1.1)

print(H.compute_average_energy_fast(N, 1))
print(H.compute_average_energy(SpinConfig("+" * N), 1))

