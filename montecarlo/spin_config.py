class SpinConfig:
    """
    1-D spin configuration
    """
    def __init__(self, spin_str, pbc=True):
        """
        Initialize instance

        Parameters
        ----------
        spin_str   : String
            set configuration to the bitstring
        pbc        : bool, default: true
            Should we use periodic boundary conditions?

        Returns
        -------
        """
        self.__spin_str = [-1 if x == '-' or x == '0' else 1 for x in spin_str]
        self.__index = 0
        self.__pbc = pbc

    def __repr__(self):
        return ''.join(['↓' if x == -1 else '↑' for x in self.__spin_str])
    
    def __getitem__(self, item):
        return self.__spin_str[item]
    
    def __len__(self):
        return len(self.__spin_str)
    
    def pbc(self):
        """
        Checks if periodic boundary conditions are in use
        
        Parameters
        ----------
        
        Returns
        -------
        pbc : bool
            True if periodic boundary conditions are in use
        """
        return self.__pbc
    
    def magnetization(self):
        """
        Return net magnetization of current configuration 
        
        Parameters
        ----------
        
        Returns
        -------
        m : float
            magnetization
        """
        return sum(self.__spin_str)