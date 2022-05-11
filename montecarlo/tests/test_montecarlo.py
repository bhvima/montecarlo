"""
Unit and regression test for the montecarlo package.
"""

import sys

import pytest
import random
from montecarlo import *

def test_montecarlo_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "montecarlo" in sys.modules

@pytest.fixture
def hamiltonian():
    """Returns an Ising Hamiltonian with k=-1, J=-2 and mu=1.1"""
    return IsingHamiltonian(k=-1, J=-2, mu=1.1)

@pytest.fixture
def spin_config_sm():
    """Returns an Spin Configuration with only 2 sites"""
    return SpinConfig('-+')

@pytest.fixture
def spin_config_without_pbc_sm():
    """Returns an Spin Configuration without Periodic Boundary Conditions with only 2 sites"""
    return SpinConfig('-+', False)

@pytest.fixture
def spin_config():
    """Returns an Spin Configuration"""
    return SpinConfig('-+--+')

@pytest.fixture
def spin_config_without_pbc():
    """Returns an Spin Configuration without Periodic Boundary Conditions"""
    return SpinConfig('-+--+', False)

def test_SpinConfigInit():
    """Tests the different methods to initialize the SpinConfig class"""
    assert list(SpinConfig('-+--++-+')) == list(SpinConfig('01001101', False))

def test_SpinConfig(spin_config, spin_config_without_pbc):
    """Tests the SpinConfig class"""
    assert list(spin_config) == [-1, 1, -1, -1, 1]
    assert str(spin_config) == '↓↑↓↓↑'
    assert spin_config.pbc() is True
    assert spin_config_without_pbc.pbc() is False
    assert len(spin_config) == 5
    assert len(spin_config_without_pbc) == 5
    assert spin_config.magnetization() == -1

def test_IsingHamiltonian(hamiltonian, spin_config_sm, spin_config_without_pbc_sm):
    """Tests the IsingHamiltonian class"""
    assert hamiltonian.energy(spin_config_sm) == 4.0
    assert hamiltonian.energy(spin_config_without_pbc_sm) == 2.0

    assert hamiltonian.compute_average_energy(spin_config_sm, 1) == -6.145889782051417
    assert hamiltonian.compute_average_magnetization(spin_config_sm, 1) == 1.9513429553841268
    assert hamiltonian.compute_heat_capacity(spin_config_sm, 1) == 0.23950137558927764
    assert hamiltonian.compute_magnetic_susceptibility(spin_config_sm, 1) == 0.19196693603718007

    assert hamiltonian.compute_average_energy_fast(2, 1) == -6.145889782051416
    assert hamiltonian.compute_average_magnetization_fast(2, 1) == 1.9513429553841268
    assert hamiltonian.compute_heat_capacity_fast(2, 1) == 0.23950137558929185
    assert hamiltonian.compute_magnetic_susceptibility_fast(2, 1) == 0.19196693603718007

    assert hamiltonian.compute_average_energy_fast(10, 1) == -30.99562217946881
    assert hamiltonian.compute_average_magnetization_fast(10, 1) == 9.999059904921909
    assert hamiltonian.compute_heat_capacity_fast(10, 1) == 0.046131068738759495
    assert hamiltonian.compute_magnetic_susceptibility_fast(10, 1) == 0.0023485385990369423

    assert hamiltonian.compute_average_energy_fast(11, 1) == -34.095184403383776
    assert hamiltonian.compute_average_magnetization_fast(11, 1) == 10.998965900847892
    assert hamiltonian.compute_heat_capacity_fast(11, 1) == 0.05074404624292583
    assert hamiltonian.compute_magnetic_susceptibility_fast(11, 1) == 0.0025832851479918872

    assert hamiltonian.compute_average_energy(spin_config_without_pbc_sm, 1) == -4.122087557865234
    assert hamiltonian.compute_average_magnetization(spin_config_without_pbc_sm, 1) == 1.943692764174632
    assert hamiltonian.compute_heat_capacity(spin_config_without_pbc_sm, 1) == 0.38131320356151477
    assert hamiltonian.compute_magnetic_susceptibility(spin_config_without_pbc_sm, 1) == 0.2060839557683165

    random.seed(3)
    assert hamiltonian.compute_average_energy_metropolis(10, 1) == -30.705599999999983
    assert hamiltonian.compute_average_magnetization_metropolis(10, 1) == 9.945999999999993
    assert hamiltonian.compute_heat_capacity_metropolis(10, 1) == 946.9159200000003
    assert hamiltonian.compute_magnetic_susceptibility_metropolis(10, 1) == 99.316
