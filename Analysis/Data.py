#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 15:17:25 2020

@author: alberto
"""

from scipy.constants import pi, c
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate



data = np.loadtxt("OSA_8.txt")

# %% Section 1 
plt.figure()
plt.plot(data[:,0], data[:,1])
plt.grid()
plt.show()

power = np.transpose(10**(data[:,1]/10)*10**-3)

plt.figure()
plt.plot(data[:,0], power)
plt.grid()
plt.show()

#%%
plt.figure()
n, bins, patches = plt.hist(power)
plt.show()

# %% Section 2

ps = np.fft.fftshift(np.fft.fft(power))
nu = c/(data[:,0]*10**-9)

plt.figure()
plt.plot(nu, np.abs(ps))
plt.show()


example = np.convolve(data[:,1], data[:,1], 'same')
plt.plot(data[:,0], example)
plt.show()



# %% Section 3 

P = integrate.quad(power, np.min(power), np.max(power))