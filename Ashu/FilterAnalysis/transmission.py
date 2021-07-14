import numpy as np
import math
c = 3e8 # speed of light in m/s

class ABCD:
    """
    Obtain S-parameters from the ABCd matrix of a 2-port network
    """
    def __init__(self, abcd, Z_l=100):
        """Load an ABCD matrix for analysis. Defnintions adopted from Pozar (2011)

        Args:
            abcd (np.ndarray): An ABCD matrix with dimension (2,2,n)
            Z_l (int, optional): The terminal impedances. It is assumed impedances 
                                 are matched at both ends. Defaults to 100.
        """
        self.abcd = abcd
        self.Z_l = Z_l
        
    def S11(self):
        A, B, C, D = (
            self.abcd[0,0], 
            self.abcd[0,1],
            self.abcd[1,0],
            self.abcd[1,1]
        )
        Z_l = self.Z_l
        num = A + B/Z_l - C*Z_l - D
        den = A + B/Z_l + C*Z_l + D
        return num/den
    
    def S12(self):
        A, B, C, D = (
            self.abcd[0,0], 
            self.abcd[0,1],
            self.abcd[1,0],
            self.abcd[1,1]
        )
        Z_l = self.Z_l
        num = 2* (A*D -B*C)
        den = A + B/Z_l + C*Z_l + D
        return num/den
    
    def S21(self):
        A, B, C, D = (
            self.abcd[0,0], 
            self.abcd[0,1],
            self.abcd[1,0],
            self.abcd[1,1]
        )
        Z_l = self.Z_l
        num = 2
        den = A + B/Z_l + C*Z_l + D
        return num/den
    
    def S22(self):
        A, B, C, D = (
            self.abcd[0,0], 
            self.abcd[0,1],
            self.abcd[1,0],
            self.abcd[1,1]
        )
        Z_l = self.Z_l
        num = -A + B/Z_l - C*Z_l + D
        den = A + B/Z_l + C*Z_l + D
        return num/den  

    
def abcd_series(Z0, B, l):
    """Product an ABCD matrix for a series transmission line

    Args:
        Z0 (float): Impedance of transmission line
        B (np.ndarray): 1-D array of wavenumbers for frequency range of length n (B = 2pi*f/c) 
        l (float): Physical length of trnasmission line

    Returns:
        np.ndarray: ABCD matrix for series trnsmission line with dims (2,2,n)
    """
    bl = B*l
    sin = np.sin(bl)
    cos = np.cos(bl)
    return np.array(
        [
            [cos, 1j*Z0*sin],
            [1j*1/Z0*sin,cos]
        ]
    )


def cascade_series(n, low, high, low_first=True, Z_l=50):
    """Cascade ABCD matrices for an alternating series of n inductors
       Note: use same array of wavenumbers for `low` and `high`       
    
    Args:
        n (int): Number of inductors
        low (np.ndarray): result from func:abcd_series. low impedance resonator.
        high (np.ndarray): result from func:abcd_series. high impedance resonator
        low_first (bool, optional): start with low impedance inductor. Defaults to True.
        Z_l (float, optional): terminal impedances. Defaults to 50.

    Returns:
        [np.ndarray]: resultant ABCD matrix
    """
    current = low_first
    n_pts = np.shape(low)[2]
    abcd_mat = np.array([
        [np.ones(n_pts) ,np.zeros(n_pts)],
        [np.zeros(n_pts),np.ones(n_pts)]
    ])
    low = np.moveaxis(low, 2, 0)
    high = np.moveaxis(high, 2, 0)
    abcd_mat = np.moveaxis(abcd_mat, 2, 0)
    for _ in range(n):
        if current:
            abcd_mat = np.matmul(abcd_mat, low)
        else:
            abcd_mat = np.matmul(abcd_mat, high)
        # switch
        current= not current
    abcd_mat = np.moveaxis(abcd_mat, 0, 2)
    return  ABCD(abcd_mat, Z_l = Z_l)


def abcd_shunt(Z0, B, l):
    """Product an ABCD matrix for a shunted transmission line

    Args:
        Z0 (float): Impedance of transmission line
        B (np.ndarray): 1-D array of wavenumbers for frequency range of length n (B = 2pi*f/c) 
        l (float): Physical length of trnasmission line

    Returns:
        np.ndarray: ABCD matrix for series trnsmission line with dims (2,2,n)
    """
    bl = B*l
    tan = np.tan(bl)
    ones = np.ones_like(tan)
    zeros = np.zeros_like(tan)
    return np.array(
        [
            [ones, zeros],
            [1j*1/Z0*tan,ones]
        ]
    )


def cascade_stub(n, shunt, resonator, Z_l=50):
    """ Cascade ABCD matrices for shunted transmission line    
    
    --------|Inductor|-------        --------
        |               |                |
     |Shunt|         |Shunt|   ....   |Shunt|
        |               |                | 
    -------------------------        --------
    Note: use same array of wavenumbers for `shunt` and `resonator`

    Args:
        n (int): Number of inductors
        shunt (np.ndarray): result of func:abcd_shunt
        resonator ([type]): result of func:abcd_series
        Z_l (float, optional): terminal impedances. Defaults to 50.

    Returns:
        [np.ndarray]: resultant ABCD matrix
    """
    n_pts = np.shape(shunt)[2]
    abcd_mat = np.array([
        [np.ones(n_pts) ,np.zeros(n_pts)],
        [np.zeros(n_pts),np.ones(n_pts)]
    ])
    shunt = np.moveaxis(shunt, 2, 0)
    resonator = np.moveaxis(resonator, 2, 0)
    abcd_mat = np.moveaxis(abcd_mat, 2, 0)
    abcd_mat = np.matmul(abcd_mat, shunt)
    for _ in range(n):
        abcd_mat = np.matmul(abcd_mat, resonator)
        abcd_mat = np.matmul(abcd_mat, shunt)
    abcd_mat = np.moveaxis(abcd_mat, 0, 2)
    return  ABCD(abcd_mat, Z_l = Z_l)