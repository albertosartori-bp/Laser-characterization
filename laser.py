import numpy as np
import scipy.integrate as integ
import scipy.optimize as opt
import matplotlib.pyplot as plt


NUM_ARRAY = [31, 26, 25, 27, 25, 14]
START = 0
DATASET = 2
NUM = NUM_ARRAY[DATASET - 1]
name = 'Measure/T_{}_{}.txt'


a = np.loadtxt(name.format(DATASET, START), skiprows=4, unpack=True)  # NOTE: for different datasets the length could be different
current = np.empty(NUM)
temperature = np.empty(NUM)
wlength = np.empty((NUM, a.shape[1]))
dbm_power = np.empty((NUM, a.shape[1])) # NOTE: every "-10dbm" is an order of magnitude less
power = np.empty((NUM, a.shape[1]))
max_back = np.empty((NUM, a.shape[1]), dtype=bool)


#%%    READ DATA, CONVERSION FROM dBm TO W, FIND MAXIMA
for i in range(NUM):
    
    # read data
    temperature[i], current[i]  = np.loadtxt(name.format(DATASET, i+START), skiprows=2, max_rows=2, unpack=True)
    wlength[i], dbm_power[i] = np.loadtxt(name.format(DATASET, i+START), skiprows=4, unpack=True)
    
    # conversion
    power[i] = 10**(dbm_power[i]/10)
    
    # finds maxima
    max_back[i] = np.r_[False, power[i,1:] > power[i,:-1]] & np.r_[power[i,:-1] > power[i,1:], False]
max_wlength = [wlength[i, np.argmax(dbm_power[i])] for i in range(NUM)]
max_p_dbm = [np.max(dbm_power[i, max_back[i]]) for i in range(NUM)]
max_p = [np.max(power[i, max_back[i]]) for i in range(NUM)]
tot_p = np.sum(power, axis=1)


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


#%% PLOT ALL SPECTRA
plt.figure()
for i in range(NUM):
    plt.plot(wlength[i], dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()






#%%  TOTAL POWER VS CURRENT STUDY

CURRENT_LIMITS_BELOW = np.array([ [0, 4],
                         [0, 7],
                         [0, 8],
                         [0, 10],
                         [0, 13],
                         [0, 0] ])

CURRENT_LIMITS_ABOVE = np.array([ [10.5, 30],
                         [11, 30],
                         [12, 30],
                         [13.5, 30],
                         [16, 30],
                         [12, 30] ])

retta = lambda x,a,b: a*x + b

mask_below = (current >= CURRENT_LIMITS_BELOW[DATASET-1, 0]) & (current<= CURRENT_LIMITS_BELOW[DATASET-1, 1] )
mask_threshold = (current >= CURRENT_LIMITS_BELOW[DATASET-1, 1]) & (current<= CURRENT_LIMITS_ABOVE[DATASET-1, 0] )
mask_above = (current >= CURRENT_LIMITS_ABOVE[DATASET-1, 0]) & (current<= CURRENT_LIMITS_ABOVE[DATASET-1, 1] )

# Linear fit below and above the threshold
param_below, covm_below = opt.curve_fit(retta, current[mask_below], tot_p[mask_below])
param_above, covm_above = opt.curve_fit(retta, current[mask_above], tot_p[mask_above])


###     FIRST METHOD: LINEAR FIT
i_th1 = -param_above[1]/param_above[0]

###     SECOND METHOD: TWO-SEGMENT FIT
i_th2 = - (param_above[1]-param_below[1])/(param_above[0]-param_below[0])

###     THIRD METHOD: FIRST DERIVATIVE
first_deriv = (tot_p[1:] - tot_p[:-1])/(current[1:] - current[:-1])
first_deriv_current = (current[1:] + current[:-1])/2

###     FOURTH METHOD: SECOND DERIVATIVE
second_deriv = (first_deriv[1:] - first_deriv[:-1])/(first_deriv_current[1:] - first_deriv_current[:-1])
second_deriv_current = (first_deriv_current[1:] + first_deriv_current[:-1])/2
i_th4 = second_deriv_current[np.argmax(second_deriv)]





#%% POWER VS CURRENT PLOT
fig, axs = plt.subplots(2,2)

axs[0,0].plot(current, tot_p, '.-')
axs[0,0].plot(np.arange(8,30), param_above[0]*np.arange(8,30) + param_above[1], '.-')
axs[0,0].grid(True)

axs[0,1].plot(current, tot_p, '.-')
axs[0,1].plot(np.arange(8,30), param_above[0]*np.arange(8,30) + param_above[1], '.-')
axs[0,1].plot(np.arange(0,15), param_below[0]*np.arange(0,15) + param_below[1], '.-')
axs[0,1].grid(True)

axs[1,0].plot(first_deriv_current, first_deriv, '.-')
axs[1,0].grid(True)

axs[1,1].plot(second_deriv_current, second_deriv, '.-')
axs[1,1].grid(True)

plt.xlabel('Current [mA]')
plt.ylabel('Optical Power [mW]')
plt.show()








#%%     STATISTICS AND SHAPE OF SPECTRUM
# Which is the wavelength of the max peak as a function of the current
plt.figure()
plt.plot(current[3:], max_wlength[3:], '.-') #skip first 3 points that are under threshold
plt.grid(True)
plt.xlabel("Current [mA]")
plt.ylabel("Wavelength [nm]")

# Count how many longitudinal modes are present. tol is the tolerance with which I select a peak
# tol = -10 means I count all the peaks that are between 0 and -10 dbm after normalization with the highest peak
# NOTE: THERE CAN BE FALSE COUNTS! A volte ci sono 2 puntine su un picco
tol = -20
norm_dbm_power = np.transpose(np.transpose(dbm_power) - max_p_dbm)
peak_count = [np.sum(norm_dbm_power[i, max_back[i]]>tol) for i in range(NUM)]
fig,ax = plt.subplots()
ax.plot(current[13:], peak_count[13:], '.-') #skip first points that are under threshold
ax.grid(True)
ax.set_xlabel("Current [mA]")
ax.set_ylabel("Longitudinal mode counts")
ax2 = ax.twinx()
ax2.plot(current[13:], tot_p[13:], '.r-') #skip first points that are under threshold
ax2.set_ylabel("Total power") #(or the one of the highest peak) non vedo nessun pattern








# # %% power/temperature plot
# plt.figure()
# plt.plot(temperature, wlength[0,np.argmax(power, axis=1)], '.-')

# plt.xlabel('Temperature [°C]')
# plt.ylabel('Peak wavelength [nm]')
# plt.grid()
# plt.show()


#%% SINGLE PLOT
roba = 1
plt.figure()
plt.plot(wlength[-roba], norm_dbm_power[-roba])
plt.grid()
plt.show()


# # %% Convolve gaussian psf of photodiode with 
# def gaussian(x, mu, sig):
#     return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# x = np.linspace(-1, 1, 200)
# y = gaussian(x, 0, 0.01)




# wlength, dbm_power = np.loadtxt('OSA_8.txt'.format(i+9), skiprows=4, unpack=True)
# power= 10**(dbm_power/10)



# h = np.convolve(power[6031:6231],y, 'same')

# plt.figure()
# plt.plot(wlength[6031:6231], h)
# plt.plot(wlength[6031:6231], power[6031:6231])
# plt.show()