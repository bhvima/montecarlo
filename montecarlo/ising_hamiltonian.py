import numpy as np
import random
from math import comb
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

    @staticmethod
    def N(n, r, k):
        return (comb(n - r - 1, k - 1) * comb(r - 1, k - 1) * n)//k

    def E(self, n, r, k):
        return (-self.J * (n - 4 * k) + self.mu * (n - 2 * r))/self.k

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
        return num/denum

    def compute_average_energy_fast(self, n, T, power=1):
        """ Compute Average Energy exactly (fast) (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        
        Returns
        -------
        E  : float 
            Energy
        """
        num = denum = 0
        for r in range(n + 1):
            if r == 0 or r == n:
                energy = self.E(n, r, 0)
                prob = np.exp(-energy/T)
                num += (energy ** power) * prob
                denum += prob
            for k in range(1, min(r, n - r) + 1):
                energy = self.E(n, r, k)
                freq = self.N(n, r, k)
                prob = np.exp(-energy/T)
                num += freq * (energy ** power) * prob
                denum += freq * prob
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

    def compute_average_magnetization_fast(self, n, T, power=1):
        """ Compute Average Magnetization exactly (fast) (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        
        Returns
        -------
        M  : float
            Magnetization
        """
        num = denum = 0
        for r in range(n + 1):
            if r == 0 or r == n:
                energy = self.E(n, r, 0)
                prob = np.exp(-energy/T)
                num += ((n - 2 * r) ** power) * prob
                denum += prob
            for k in range(1, min(r, n - r) + 1):
                energy = self.E(n, r, k)
                freq = self.N(n, r, k)
                prob = np.exp(-energy/T)
                num += freq * ((n - 2 * r) ** power) * prob
                denum += freq * prob
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

    def compute_heat_capacity_fast(self, n, T):
        """ Compute Heat Capacity exactly (fast) (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        
        Returns
        -------
        HC : float
            Heat Capacity
        """
        return (self.compute_average_energy_fast(n, T, power=2) - (self.compute_average_energy_fast(n, T) ** 2)) / (T ** 2)

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

    def compute_magnetic_susceptibility_fast(self, n, T):
        """ Compute Magnetic Susceptibility exactly (fast) (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        
        Returns
        -------
        MS : float
            Magnetic Susceptability
        """
        return (self.compute_average_magnetization_fast(n, T, power=2) - (self.compute_average_magnetization_fast(n, T) ** 2)) / T

    def metropolis_sweep(self, conf, T):
        """Perform a single sweep through all the sites and return updated configuration (pbc)
        Parameters
        ----------
        conf   : :class:`SpinConfig`
            input configuration 
        T      : int
            Temperature
        
        Returns
        -------
        conf  : :class:`SpinConfig`
            Returns updated config
        """
        assert(conf.pbc())
        
        N = len(conf)
        for i, v in enumerate(conf):
            delta_e = ((-2 * v * self.mu) + (-2 * v * (conf[(i - 1) % N] + conf[(i + 1) % N])))/self.k

            if delta_e <= 0 or random.random() <= np.exp(-delta_e/T):
                conf[i] = -v

        return conf

    def compute_average_energy_metropolis(self, n, T, nsweep=1000, nburn=100):
        """ Compute Average Energy approximate (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        nsweep : int, default: 1000
            monte carlo steps
        nburn  : int, default: 100
            burn steps
        
        Returns
        -------
        E  : float 
            Energy
        """
        conf_str = ["+"] * (n//2) + ["-"] * (n - n//2)
        random.shuffle(conf_str)
        conf = SpinConfig(''.join(conf_str))

        for _ in range(nburn):
            self.metropolis_sweep(conf, T=T)

        self.metropolis_sweep(conf, T=T)

        E = self.energy(conf)

        for i in range(1,nsweep):
            self.metropolis_sweep(conf, T=T)
            E = (E * i + self.energy(conf))/(i + 1)

        return E

    def compute_average_magnetization_metropolis(self, n, T, nsweep=1000, nburn=100):
        """ Compute Average Magnetization approximate (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        nsweep : int, default: 1000
            monte carlo steps
        nburn  : int, default: 100
            burn steps

        Returns
        -------
        M  : float
            Magnetization
        """
        conf_str = ["+"] * (n//2) + ["-"] * (n - n//2)
        random.shuffle(conf_str)
        conf = SpinConfig(''.join(conf_str))

        for _ in range(nburn):
            self.metropolis_sweep(conf, T=T)

        self.metropolis_sweep(conf, T=T)

        M = conf.magnetization()

        for i in range(1,nsweep):
            self.metropolis_sweep(conf, T=T)
            M = (M * i + conf.magnetization())/(i + 1)

        return M

    def compute_heat_capacity_metropolis(self, n, T, nsweep=1000, nburn=100):
        """ Compute Heat Capacity approximate (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        nsweep : int, default: 1000
            monte carlo steps
        nburn  : int, default: 100
            burn steps
        
        Returns
        -------
        HC : float
            Heat Capacity
        """
        conf_str = ["+"] * (n//2) + ["-"] * (n - n//2)
        random.shuffle(conf_str)
        conf = SpinConfig(''.join(conf_str))

        for _ in range(nburn):
            self.metropolis_sweep(conf, T=T)

        self.metropolis_sweep(conf, T=T)

        E = self.energy(conf)
        EE = E * E

        for i in range(1,nsweep):
            self.metropolis_sweep(conf, T=T)
            Ei = self.energy(conf)
            E = (E * i + Ei)/(i + 1)
            EE = (EE * i + Ei*Ei)/(i + 1)

        return (EE-E*E)/(T*T)

    def compute_magnetic_susceptibility_metropolis(self, n, T, nsweep=1000, nburn=100):
        """ Compute Magnetic Susceptibility approximate (pbc)

        Parameters
        ----------
        n      : int
            number of sites
        T      : int
            Temperature
        nsweep : int, default: 1000
            monte carlo steps
        nburn  : int, default: 100
            burn steps
        
        Returns
        -------
        MS : float
            Magnetic Susceptibility
        """
        conf_str = ["+"] * (n//2) + ["-"] * (n - n//2)
        random.shuffle(conf_str)
        conf = SpinConfig(''.join(conf_str))

        for _ in range(nburn):
            self.metropolis_sweep(conf, T=T)

        self.metropolis_sweep(conf, T=T)

        M = conf.magnetization()
        MM = M * M

        for i in range(1,nsweep):
            self.metropolis_sweep(conf, T=T)
            Mi = conf.magnetization()
            M = (M * i + Mi)/(i + 1)
            MM = (MM * i + Mi*Mi)/(i + 1)

        return (MM-M*M)/T
