Getting Started
===============

This page details how to get started with montecarlo. 


Basic Installation
------------------

To install montecarlo, run:

.. code-block:: bash
        
        git clone https://github.com/bhvima/montecarlo.git
        (cd montecarlo; pip install .)


Examples
--------

Compute the energy of a given spin configuration:

.. code-block:: python
        
        from montecarlo import IsingHamiltonian as Hamiltonian
        from montecarlo import SpinConfig

        H = Hamiltonian(k=1, J=2, mu=1.1)
        H.energy(SpinConfig("++-+---+--+"))


Plot of observables vs T
------------------------

The following is produced using the exact computation (fast version):

.. code-block:: python

	from montecarlo import IsingHamiltonian as Hamiltonian
	import matplotlib.pyplot as plt
	import numpy as np

	H = Hamiltonian(k=1, J=-2, mu=1.1)
	T = np.linspace(.1,10,100)

	AE = H.compute_average_energy_fast(30, T)
	AM = H.compute_average_magnetization_fast(30, T)
	HC = H.compute_heat_capacity_fast(30, T)
	MS = H.compute_magnetic_susceptibility_fast(30, T)

	fig = plt.figure()
	plt.plot(T,AE,'r', label='Average Energy')
	plt.plot(T,AM,'b', label='Average Magnetization')
	plt.plot(T,HC,'g', label='Heat Capacity')
	plt.plot(T,MS,'y', label='Magnetic Susceptibility')

	plt.legend()
	plt.show()

This should produce the following plot:

.. image:: ./plot1.png
  :width: 400
  :alt: Observables vs Temperature (30 sites)

The following is produced using the appoximate computation using metropolis sampling:

.. code-block:: python

	from montecarlo import IsingHamiltonian as Hamiltonian
	import matplotlib.pyplot as plt
	import random

	H = Hamiltonian(k=1, J=-2, mu=1.1)

	random.seed(3)

	AE = []
	AM = []
	HC = []
	MS = []
	T = []

	for i in range(1,100):
	    T.append(.1*i)
	    AE.append(H.compute_average_energy_metropolis(30, t[-1]))
	    AM.append(H.compute_average_magnetization_metropolis(30, t[-1]))
	    HC.append(H.compute_heat_capacity_metropolis(30, t[-1], nsweep=4000, nburn=1000))
	    MS.append(H.compute_magnetic_susceptibility_metropolis(30, t[-1], nsweep=4000, nburn=1000))

	plt.plot(T,AE,'r', label='Average Energy')
	plt.plot(T,AM,'b', label='Average Magnetization')
	plt.plot(T,HC,'g', label='Heat Capacity')
	plt.plot(T,MS,'y', label='Magnetic Susceptibility')

	plt.legend()
	plt.show()

This should produce the following plot:

.. image:: ./plot2.png
  :width: 400
  :alt: Observables vs Temperature (30 sites)