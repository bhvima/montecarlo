import numpy as np
from .spin_config import SpinConfig

class IsingHamiltonian:
    """Class for Ising Hamiltonian

        .. math::
            H = -J\\sum_{\\left<ij\\right>} \\sigma_i\\sigma_j + \\mu\\sum_i\\sigma_i

    """

    def __init__(self, k, J=1.0, mu=0.0):
        """ Constructor

        Parameters
        ----------
        k: float, Boltzmann constant
        J: float
            Strength of coupling
        mu: float
            Chemical potential
        """
        self.k = k
        self.J = J
        self.mu = mu

    def energy(self, config):
        """Compute energy of configuration, `config`

            .. math::
                E = \\left<\\hat{H}\\right>

        Parameters
        ----------
        config   : :class:`SpinConfig`
            input configuration

        Returns
        -------
        energy  : float
            Energy of the input configuration
        """
        if config.pbc():
            return (-self.J * sum(config[i - 1] * config[i] for i in range(len(config))) + self.mu * sum(config))/self.k
        return (-self.J * sum(config[i] * config[i + 1] for i in range(len(config) - 1)) + self.mu * sum(config))/self.k


    def __probability(self, config, T):
        return np.exp((-1/T) * self.energy(config))

    def compute_average_energy(self, conf, T, power=1):
        """ Compute Average Energy exactly

        Parameters
        ----------
        conf   : :class:`SpinConfig`
            input configuration
        T      : int
            Temperature
        
        Returns
        -------
        E  : float 
            Energy
        """
        num = denum = 0
        for x in range(2**len(conf)):
            s = SpinConfig(bin(x)[2:].zfill(len(conf)), conf.pbc())
            p = self.__probability(s, T)
            num += (self.energy(s) ** power) * p
            denum += p
        print(denum)
        return num/denum

    def compute_average_magnetization(self, conf, T, power=1):
        """ Compute Average Magnetization exactly

        Parameters
        ----------
        conf   : :class:`SpinConfig`
            input configuration 
        T      : int
            Temperature
        
        Returns
        -------
        M  : float
            Magnetization
        """
        num = denum = 0
        for x in range(2**len(conf)):
            s = SpinConfig(bin(x)[2:].zfill(len(conf)), conf.pbc())
            p = self.__probability(s, T)
            num += (s.magnetization() ** power) * p
            denum += p
        return num/denum

    def compute_heat_capacity(self, conf, T):
        """ Compute Heat Capacity exactly

        Parameters
        ----------
        conf   : :class:`SpinConfig`
            input configuration 
        T      : int
            Temperature
        
        Returns
        -------
        HC : float
            Heat Capacity
        """
        return (self.compute_average_energy(conf, T, power=2) - (self.compute_average_energy(conf, T) ** 2)) / (T ** 2)

    def compute_magnetic_susceptibility(self, conf, T):
        """ Compute Magnetic Susceptibility exactly

        Parameters
        ----------
        conf   : :class:`SpinConfig`
            input configuration 
        T      : int
            Temperature
        
        Returns
        -------
        MS : float
            Magnetic Susceptability
        """
        return (self.compute_average_magnetization(conf, T, power=2) - (self.compute_average_magnetization(conf, T) ** 2)) / T
