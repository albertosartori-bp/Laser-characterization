import numpy as np
import scipy.integrate as integ
import scipy.optimize as opt
import matplotlib.pyplot as plt

NUM = 17
START = 0

a = np.loadtxt("Measure/I_1_0.txt", skiprows=4, unpack=True)  # NOTE: for different datasets the length could be different
current = np.empty(NUM)
temperature = np.empty(NUM)
wlength = np.empty((NUM, a.shape[1]))
dbm_power = np.empty((NUM, a.shape[1])) # NOTE: every "-10dbm" is an order of magnitude less
power = np.empty((NUM, a.shape[1]))
dpower = np.empty((NUM, a.shape[1], 2))
max_back = np.empty((NUM, a.shape[1]), dtype=bool)
dtot_p = np.empty((NUM,2))

#%%    READ DATA, CONVERSION FROM dBm TO W, FIND MAXIMA
for i in range(NUM):
    
    # read data
    temperature[i], current[i]  = np.loadtxt("Measure/I_1_{}.txt".format(i+START), skiprows=2, max_rows=2, unpack=True)
    wlength[i], dbm_power[i] = np.loadtxt("Measure/I_1_{}.txt".format(i+START), skiprows=4, unpack=True)
    
    # conversion
    power[i] = 10**(dbm_power[i]/10)/3.8
    dpower[i,:,0] = power[i] * 10**(-0.5/10) - power[i]
    dpower[i,:,1] = power[i] * 10**(0.5/10) - power[i]
    
    # finds maxima
    max_back[i] = np.r_[False, power[i,1:] > power[i,:-1]] & np.r_[power[i,:-1] > power[i,1:], False]
max_wlength = [wlength[i, np.argmax(dbm_power[i])] for i in range(NUM)]
max_p_dbm = [np.max(dbm_power[i, max_back[i]]) for i in range(NUM)]
max_p = [np.max(power[i, max_back[i]]) for i in range(NUM)]

tot_p = np.sum(power, axis=1)
dtot_p[:, 0] = np.sqrt(np.sum(np.squeeze(dpower[:,:,0])**2, axis=1))
dtot_p[:, 1] = np.sqrt(np.sum(np.squeeze(dpower[:,:,1])**2, axis=1))
#print(dtot_p/np.transpose([tot_p, tot_p]))

# #%% PLOT TEMPERATURE AND CURRENT
# fig, axs = plt.subplots(2)
# axs[0].plot(temperature, '.-g', label="Temperature")
# axs[0].grid(True)
# axs[0].set(ylabel = "Temperature [°C]")
# axs[1].plot(current, ".-", label="Current")
# plt.xlabel('Different acquisitions')
# plt.ylabel("Current [mA]")
# axs[1].grid(True)
# plt.show()


# #%% PLOT ALL SPECTRA
# plt.figure()
# for i in range(NUM):
#     plt.plot(wlength[i], dbm_power[i], label='{}'.format(i+1))
#     #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

# plt.legend(loc='upper right')
# plt.xlabel('wavelength [nm]')
# plt.ylabel('power [dBm]')
# plt.grid(True)
# plt.show()


#%% PLOT ONE SPECTRUM
plt.figure()
i = 3
plt.plot(wlength[i], dbm_power[i], '.-', label='{}'.format(i+1))
plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], '.', label='{}'.format(i+1))

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()



# %% power/temperature plot

#find the first N maxima (Nmax)
Nmax = 5
N = 15
idx = np.empty((NUM, N))
y_temp = np.empty((NUM, N))
y = np.empty((NUM, Nmax))
for i in range(NUM):
    idx[i] = (-power[i, max_back[i]]).argsort()[:N]
idx = idx.astype(int)
for i in range(NUM):
    y_temp[i] = wlength[i, max_back[i]][idx[i]]
    
    
# clean y: some peaks have multiple maxima
for i in range(NUM):
    j=0
    k = 0
    while j<5:
        if np.all( abs(y_temp[i,:k] - y_temp[i,k]) > 0.7):
            y[i,j] = y_temp[i, k]
            j += 1
        k += 1

# fit
retta = lambda x,a,b: a*x + b
fit_param, fit_covm = opt.curve_fit(retta, temperature[:8], y[:8, 0], sigma = 0.02*np.ones(8), absolute_sigma = True)
yerr = np.sqrt( (0.02*np.ones(8))**2 + (fit_param[0]*0.1)**2 ) # RIVEDI ERRORE SU TEMPERATURA
fit_param, fit_covm = opt.curve_fit(retta, temperature[:8], y[:8, 0], sigma = yerr, absolute_sigma = True)

# plot
plt.figure()
for i in range(Nmax):
    plt.errorbar(temperature, y[:,i], yerr = 0.02, fmt='.-', label="Maximum number {}".format(i+1))
for i in range(10):
    plt.plot(np.arange(10,50,1), fit_param[1]-8.8+1.1*i + fit_param[0]*np.arange(10,50,1), '-k', alpha=0.5)

plt.xlabel('Temperature [°C]')
plt.ylabel('Peak wavelength [nm]')
plt.grid()
plt.legend()
plt.show()

print("The coefficient is: ", fit_param[0], "+-", np.sqrt(fit_covm[0,0]), "nm/K")