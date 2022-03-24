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
        
        from montecarlo import IsingHamiltonian1D as hamiltonian
        from montecarlo import SpinConfig

        H = hamiltonian(k=1, J=2, mu=1.1)
        H.energy(SpinConfig("++-+---+--+"))
