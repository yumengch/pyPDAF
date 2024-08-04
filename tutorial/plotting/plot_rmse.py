#!/usr/bin/env python3

# A script to plot the 2D field from the online_2D_serialmodel tutorial.
# Requires Python 3, Matplotlib and Numpy.

# Usage: ./plot_rmse.py <experiment directory name> 

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

def rmse(filename1, filename2):
    field1 = np.loadtxt(filename1, delimiter=';')
    field2 = np.loadtxt(filename2)
    field1 = field1.reshape(18,36)
    field2 = field2.reshape(18,36)
    rmse = 0;
    for i in range(16):
       for j in range(36):
	       rmse = rmse + (field1[i,j]-field2[i,j])**2
    rmse = np.sqrt(1/(18*36)*(rmse));

    return rmse

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('experiment')
    args = parser.parse_args()

    nsteps = ["2","4","6","8","10","12","14","16","18"]
    
    rmse_arr = np.full((len(nsteps)), np.nan)
    for j in range(len(nsteps)):
        rmse_arr[j] = rmse(f'{args.experiment}/state_step{nsteps[j]}_ana.txt', f'inputs_online/true_step{nsteps[j]}.txt')

    plt.title(f'Experiment: {args.experiment}', fontsize=18)
    plt.xlabel('analysis time step',fontsize=16)
    plt.ylabel('RMSE', fontsize=16)

    plt.plot(nsteps,rmse_arr)
    plt.show()
