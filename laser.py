import numpy as np
import scipy.integrate as integ
import scipy.optimize as opt
import matplotlib.pyplot as plt

#%% Initialization
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
dpower = np.empty((NUM, a.shape[1], 2))
max_back = np.empty((NUM, a.shape[1]), dtype=bool)
dtot_p = np.empty((NUM,2))
n_ext = np.empty(len(NUM_ARRAY-1))

#%%    READ DATA, CONVERSION FROM dBm TO W, FIND MAXIMA
for i in range(NUM):
    
    # read data
    temperature[i], current[i]  = np.loadtxt(name.format(DATASET, i+START), skiprows=2, max_rows=2, unpack=True)
    wlength[i], dbm_power[i] = np.loadtxt(name.format(DATASET, i+START), skiprows=4, unpack=True)
    
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



#%% WAVELENGTH DEPENDENCE FROM CURRENT (DONE THE SAME WAY AS TEMPERATURE)

#find the first N maxima (Nmax)
Nmax = 1
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
    while j<Nmax:
        if np.all( abs(y_temp[i,:k] - y_temp[i,k]) > 0.7):
            y[i,j] = y_temp[i, k]
            j += 1
        k += 1

# # fit
# retta = lambda x,a,b: a*x + b
# fit_param, fit_covm = opt.curve_fit(retta, temperature[:8], y[:8, 0], sigma = 0.02*np.ones(8), absolute_sigma = True)
# yerr = np.sqrt( (0.02*np.ones(8))**2 + (fit_param[0]*0.1)**2 ) # RIVEDI ERRORE SU TEMPERATURA
# fit_param, fit_covm = opt.curve_fit(retta, temperature[:8], y[:8, 0], sigma = yerr, absolute_sigma = True)

# plot
plt.figure()
for i in range(Nmax):
    plt.errorbar(current, y[:,i], yerr = 0.02, fmt='.-', label="Maximum number {}".format(i+1))
# for i in range(10):
#     plt.plot(np.arange(10,50,1), fit_param[1]-8.8+1.1*i + fit_param[0]*np.arange(10,50,1), '-k', alpha=0.5)

plt.xlabel('current [mA]')
plt.ylabel('Peak wavelength [nm]')
plt.grid()
plt.legend()
plt.show()

# print("The coefficient is: ", fit_param[0], "+-", np.sqrt(fit_covm[0,0]), "nm/K")





#%% WAVELENGTH DEPENDENCE FROM CURRENT

#choose 10 maxima to study 
Nmax = 10
N = 55
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
    while j<Nmax:
        if np.all( abs(y_temp[i,:k] - y_temp[i,k]) > 0.7):
            y[i,j] = y_temp[i, k]
            j += 1
        k += 1

# plot
plt.figure()
for i in range(Nmax):
    plt.errorbar(current, y[:,i], yerr = 0.02, fmt='.', label="Maximum number {}".format(i+1))

plt.xlabel('current [mA]')
plt.ylabel('Peak wavelength [nm]')
plt.grid()
plt.legend()
plt.show()

# select a peak
lamb = 1545.1
y_fit = np.zeros(NUM)
for i in range(NUM):
    if not np.all(np.abs(y[i,:]-lamb)>0.7):
        y_fit[i] = y[i,np.abs(y[i,:]-lamb)<0.7]
        
    
# fit
retta = lambda x,a,b: a*x + b
fit_param, fit_covm = opt.curve_fit(retta, current[-11:], y_fit[-11:], sigma = 0.02*np.ones(11), absolute_sigma = True)
yerr = np.sqrt( (0.02*np.ones(11))**2 + (fit_param[0]*0.1)**2 ) # RIVEDI ERRORE SU CORRENTE
fit_param, fit_covm = opt.curve_fit(retta, current[-11:], y_fit[-11:], sigma = yerr, absolute_sigma = True)

plt.figure()
plt.plot(np.arange(10,30,1), fit_param[1] + fit_param[0]*np.arange(10,30,1), '-k', alpha=0.5)
plt.errorbar(current, y_fit, yerr = 0.02, fmt='.', label="Maximum number {}".format(i+1))

print("The coefficient is: ", fit_param[0], "+-", np.sqrt(fit_covm[0,0]), "nm/mA")



#%%     STATISTICS AND SHAPE OF SPECTRUM

# Count how many longitudinal modes are present. tol is the tolerance with which I select a peak
# tol = -10 means I count all the peaks that are between 0 and -10 dbm after normalization with the highest peak
# NOTE: THERE CAN BE FALSE COUNTS! A volte ci sono 2 puntine su un picco
tol = -20
norm_dbm_power = np.transpose(np.transpose(dbm_power) - max_p_dbm)
peak_count = [np.sum(norm_dbm_power[i, max_back[i]]>tol) for i in range(NUM)]
# peak_count1 è il conteggio togliendo i picchi vicini tra loro (nota: i conteggi da 1 vengono azzerati)
peak_count1 = peak_count - np.array([np.sum( abs( wlength[i,max_back[i]][norm_dbm_power[i, max_back[i]]>tol] - np.roll(wlength[i,max_back[i]][norm_dbm_power[i, max_back[i]]>tol], 1) ) < 0.6   ) for i in range(NUM)])
fig,ax = plt.subplots()
ax.plot(current[13:], peak_count[13:], '.-') #skip first points that are under threshold
ax.plot(current[13:], peak_count1[13:], '.-') #skip first points that are under threshold
ax.grid(True)
ax.set_xlabel("Current [mA]")
ax.set_ylabel("Longitudinal mode counts")
ax2 = ax.twinx()
ax2.plot(current[13:], tot_p[13:], '.r-') #skip first points that are under threshold
ax2.set_ylabel("Total power") #(or the one of the highest peak) non vedo nessun pattern




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
mask_above = (current >= CURRENT_LIMITS_ABOVE[DATASET-1, 0]) & (current<= CURRENT_LIMITS_ABOVE[DATASET-1, 1] ) &  (current != 16.92)

# Linear fit below and above the threshold
param_below, covm_below = opt.curve_fit(retta, current[mask_below], tot_p[mask_below], sigma = np.sum(dtot_p[mask_below,:], axis=1)/2, absolute_sigma = True)
param_above, covm_above = opt.curve_fit(retta, current[mask_above], tot_p[mask_above], sigma = np.sum(dtot_p[mask_above,:], axis=1)/2, absolute_sigma = True)


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


#%% External quantum efficiency
n_ext[DATASET-1] = 0.806544*param_above[0]*np.max(y)
print(param_above)

#%% POWER VS CURRENT PLOT
fig, axs = plt.subplots(2,2)

axs[0,0].fill_between(current, tot_p - dtot_p[:,0], tot_p + dtot_p[:,1], alpha=0.4)
axs[0,0].plot(current, tot_p, '.-')
axs[0,0].plot(np.arange(8,30), param_above[0]*np.arange(8,30) + param_above[1], '.-')
axs[0,0].grid(True)

axs[0,1].fill_between(current, tot_p - dtot_p[:,0], tot_p + dtot_p[:,1], alpha=0.4)
axs[0,1].plot(current, tot_p, '.-')
axs[0,1].plot(np.arange(8,30), param_above[0]*np.arange(8,30) + param_above[1], '.-')
axs[0,1].plot(np.arange(0,15), param_below[0]*np.arange(0,15) + param_below[1], '.-')
axs[0,1].grid(True)

axs[1,0].plot(current, tot_p, '.-')
axs[1,0].plot(first_deriv_current, first_deriv, '.-')
axs[1,0].grid(True)

axs[1,1].plot(current, tot_p, '.-')
axs[1,1].plot(second_deriv_current, second_deriv, '.-')
axs[1,1].grid(True)

plt.xlabel('Current [mA]')
plt.ylabel('Optical Power [mW]')
plt.show()




#%% Characteristic temperature : threshold current as function of temperature.

i_th = [8.640025086872926, 10.515240037711957,  11.67158638894386, 13.063864899267069, 15.499155841615838]
T = np.empty(5)
for j in range(5):
        T[j] = np.loadtxt(name.format(j+1, 0), skiprows=2, max_rows=1, unpack=True)
print(T)

plt.figure()
plt.plot(T, np.log(i_th), '.r')
plt.xlabel('Temperature [°C]')
plt.ylabel('$ln(J_{th})$')
plt.grid(True)
plt.show()

param, covm = opt.curve_fit(retta, T, np.log(i_th))
x = np.arange(15,50, 0.1)

plt.plot(x, param[0]*x+param[1], '-b')
plt.legend('Experimental data', 'Fit')

#%% External quantum efficiency


plt.figure()
plt.plot(T, n_ext/4, '.r')
plt.grid(True)
plt.show()

#%%
