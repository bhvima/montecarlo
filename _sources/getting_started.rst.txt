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
