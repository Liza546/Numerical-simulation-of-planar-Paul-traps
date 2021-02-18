# put your functions here
# calculate 1D trap depth
import matplotlib.pyplot  as plt, numpy as np, scipy.constants as ct
from scipy.signal import argrelextrema

def trap_depth(s,axis="z", x=np.linspace(-2, 2, 1000), min_pos=[0, 0, 1], q=1, m=1, u=1, o=1, l=1, mode="J",
               plotting=True):
    """Function calculates the 1D trap depth and the trap size for a given axis.

    Input parameters:
    axis: "x", "y" or "z", selection of the axis along which you want to compute
    the trapping depth
    q: elementary charge
    m: mass of the particle
    u: rf-voltage amplitude
    o: rf-voltage frequency
    l: length scale
    mode: "J" for selection of calculation in Joules or "K" for selection of
    calculation in Kelvins
    plotting: "True" for plotting the graph of potential vs. axis, otherwise
    "False" without the plot
    Output parameters:
    trap depth: value of the 1D trapping depth in eV
    """
    xmin, ymin, zmin = min_pos
    p = []

    # calculate the potential at this set of points
    for n in range(len(x)):
        if axis == 'z':
            p.append(s.potential([xmin, ymin, x[n]])[0] / u ** 2)
            x_label = 'z/l'
        elif axis == 'y':
            p.append(s.potential([xmin, x[n], zmin])[0] / u ** 2)
            x_label = 'y/l'
        elif axis == 'x':
            p.append(s.potential([x[n], ymin, zmin])[0] / u ** 2)
            x_label = 'x/l'

            # create arrays from potential and z-data
    psi = np.array(p)
    # calculate extremes of psi in z
    x_max = argrelextrema(psi, np.greater, order=2)[0]  # maximum
    x_min = argrelextrema(psi, np.less, order=2)  # minimum

    if (not x_min) | (not x_max[0]):
        print('Warning: minimum/maximum point not found.')
    else:
        # calculate the trap depth for the chosen axis
        D = (psi[x_max[0]] - psi[x_min[0]]) * (q ** 2 * u ** 2) / (4 * m * o ** 2 * l ** 2)  # [J]
        delta_x = (x[x_max[0]] - x[x_min[0]]) * l  # [m]
        xmin = x[x_min[0]] * l  # [m]

    if plotting == True:
        plt.figure()
        plt.plot(x, psi)
        plt.xlabel(x_label)
        plt.ylabel('Potential (units of (Q**2*Vrf**2)/(4*m*o**2*l**2))')
        #
    if mode == 'J':
        return D, delta_x, xmin
    elif mode == 'K':
        return D / ct.k, delta_x, xmin
