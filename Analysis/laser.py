import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('qt5agg')



# %% 
x = np.linspace(2, 10, 10)
y = np.sinc(x)

plt(x,y)
# %%
NUM = 9

a = np.loadtxt('Reso_0.txt', skiprows=4, unpack=True)

current = np.empty(NUM)
temperature = np.empty(NUM)
wlength = np.empty((NUM, a.shape[1]))
dbm_power = np.empty((NUM, a.shape[1]))
power = np.empty((NUM, a.shape[1]))
max_back = np.empty((NUM, a.shape[1]), dtype=bool)


for i in range(NUM):
    
    current, temperature = np.loadtxt('Reso_{}.txt'.format(i), skiprows=2, 
                                      max_rows=2, unpack=True)
    
    wlength[i], dbm_power[i] = np.loadtxt('Reso_{}.txt'.format(i), skiprows=4, unpack=True)

    power[i] = 10**(dbm_power[i]/10)

    max_back[i] = np.r_[False, power[i,1:] > power[i,:-1]] & np.r_[power[i,:-1] > power[i,1:], False]


# %% plot all spectra
plt.figure()
for i in range(NUM):
    plt.plot(wlength[i], dbm_power[i], label='{}'.format(i+1))
    plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend()
plt.grid(True)
plt.show()

# %%  power

tot_p = integ.simps(10**(power/10), wlength)
max_p = np.max(10**(power[:,max_back[1]]/10), axis=1)

# %% Convolve gaussian psf of photodiode with 
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x = np.linspace(-1, 1, 200)
y = gaussian(x, 0, 0.01)




wlength, dbm_power = np.loadtxt('OSA_8.txt'.format(i+9), skiprows=4, unpack=True)
power= 10**(dbm_power/10)



h = np.convolve(power[6031:6231],y, 'same')

plt.figure()
plt.plot(wlength[6031:6231], h)
plt.plot(wlength[6031:6231], power[6031:6231])
plt.show()
