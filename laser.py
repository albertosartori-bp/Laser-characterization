import numpy as np
import scipy.integrate as integ
import scipy.optimize as opt
import matplotlib.pyplot as plt

i_th = np.empty([5,4])
di_th = np.empty([5,4])

e_qe = np.empty(5)
d_eqe = np.empty(5)

#%% Initialization
NUM_ARRAY = [31, 26, 25, 27, 25, 14]
START = 0
DATASET = 5
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



n_ext = np.empty(len(NUM_ARRAY)-1)

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
max_wlength = np.array([wlength[i, np.argmax(dbm_power[i])] for i in range(NUM)])
max_p_dbm = np.array([np.max(dbm_power[i, max_back[i]]) for i in range(NUM)])
max_p = np.array([np.max(power[i, max_back[i]]) for i in range(NUM)])
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
        y_fit= y[i,np.abs(y[i,:]-lamb)<0.7]
        
    
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
retta = lambda x, a, b: a*x + b


ith = [10, 10.52,  11.67, 13.06, 15.6]
mask = np.array([current[i] > ith[DATASET-1] for i in range(NUM)]) #take only the points above the threshold
mask_below = np.array([current[i] < ith[DATASET-1] for i in range(NUM)])


#exclude specific data points
if DATASET == 1:
    mask[[22,25,29]] = False  #16 e 17 can be excluded but doesn't change
elif DATASET == 2:
    mask[17] = False
elif DATASET == 3:
    mask[[18,23]] = False
elif DATASET == 4:
    mask[[25]] = False #23 e 24 can be excluded but doesn't change
elif DATASET == 5:
    mask[[16,17]] = False


param, covm = opt.curve_fit(retta, current[mask], tot_p[mask], sigma= tot_p[mask]*0.05)
print("a: {} +/- {} mW/mA".format(param[0], np.sqrt(covm[0,0])))
print("b: {} +/- {} mW".format(param[1], np.sqrt(covm[1,1])))

plt.figure()
plt.errorbar(current, tot_p, tot_p*0.05, fmt='.', capsize=2, label='Experimental data')
#plt.fill_between(current, tot_p - dtot_p[:,0], tot_p + dtot_p[:,1], alpha=0.4) #with error from datasheet
#plt.fill_between(current, tot_p*0.95, tot_p*1.05, alpha=0.4) #with error 3% (+-1.5?)
plt.plot(current[mask], retta(current[mask], *param,), label='Linear fit')
plt.xlabel('Current [mA]')
plt.ylabel('Optical power [mW]')
plt.legend()
plt.grid(True)
plt.show()

eqe = 1.6022e-19 / 6.6262e-34 / 2.99e8 *param[0]* np.mean(max_wlength[mask]) *1e-9
print('External Qefficiency: {:%} +/- {:%}'.format(eqe, eqe*np.sqrt(covm[0,0]/param[0]**2 + np.std(max_wlength[mask])**2/np.mean(max_wlength[mask])**2)))#PROPAGA ERRORE DA LAMBDA

#%% quantum efficiency in function of temperature
#slope efficiency
#0.094(2)  0.089(1)  0.082(2)  0.079(1)  0.074(1)
#external quantum efficiency
#0.117(3) 0.111(1) 0.103(2) 0.099(1) 0.093(2)
t_t =np.array([16, 23, 30, 37, 48])
t_eqe = np.array([0.117, 0.111, 0.103, 0.099, 0.093])
t_deqe = np.array([0.003, 0.001, 0.002, 0.001, 0.002])

t_32 = lambda x,a,b: b+a*x**(-5/2)

t_param, t_covm = opt.curve_fit(t_32, t_t, t_eqe, sigma=t_deqe)

plt.figure()
plt.errorbar(t_t, t_eqe, yerr=t_deqe, fmt='.', capsize=2)
plt.plot(t_t, t_32(t_t, *t_param,))
plt.xlabel('Temperature [Celsius]')
plt.ylabel('External quantum efficiency')
plt.grid(True)
plt.show()

#%% THRESHOLD CURRENT STUDY
CURRENT_LIMITS_BELOW = np.array([ [0, 4],
                         [0, 7],
                         [0, 8],
                         [0, 10],
                         [0, 13],
                         [0, 0] ])
                         

#exclude specific data points
if DATASET == 1:
    mask[[22,25,29]] = False  #16 e 17 can be excluded but doesn't change
elif DATASET == 2:
    mask[17] = False
elif DATASET == 3:
    mask[[18,23]] = False
elif DATASET == 4:
    mask[[25]] = False #23 e 24 can be excluded but doesn't change
elif DATASET == 5:
    mask[[16,17]] = False


param, covm = opt.curve_fit(retta, current[mask], tot_p[mask], sigma= tot_p[mask]*0.05)
print("a: {} +/- {} mW/mA".format(param[0], np.sqrt(covm[0,0])))
print("b: {} +/- {} mW".format(param[1], np.sqrt(covm[1,1])))

plt.figure()
plt.errorbar(current, tot_p, tot_p*0.05, fmt='.', capsize=2, label='Experimental data')
#plt.fill_between(current, tot_p - dtot_p[:,0], tot_p + dtot_p[:,1], alpha=0.4) #with error from datasheet
#plt.fill_between(current, tot_p*0.95, tot_p*1.05, alpha=0.4) #with error 3% (+-1.5?)
plt.plot(current[mask], retta(current[mask], *param,), label='Linear fit')
plt.xlabel('Current [mA]')
plt.ylabel('Optical power [mW]')
plt.legend()
plt.grid(True)
plt.show()

eqe = 1.6022e-19 / 6.6262e-34 / 2.99e8 *param[0]* np.mean(max_wlength[mask]) *1e-9
print('External Qefficiency: {:%} +/- {:%}'.format(eqe, eqe*np.sqrt(covm[0,0]/param[0]**2 + np.std(max_wlength[mask])**2/np.mean(max_wlength[mask])**2)))#PROPAGA ERRORE DA LAMBDA

#%% quantum efficiency in function of temperature
#slope efficiency
#0.094(2)  0.089(1)  0.082(2)  0.079(1)  0.074(1)
#external quantum efficiency
#0.117(3) 0.111(1) 0.103(2) 0.099(1) 0.093(2)
t_t =np.array([16, 23, 30, 37, 48])
t_eqe = np.array([0.117, 0.111, 0.103, 0.099, 0.093])
t_deqe = np.array([0.003, 0.001, 0.002, 0.001, 0.002])

t_32 = lambda x,a,b: b+a*x**(-5/2)

t_param, t_covm = opt.curve_fit(t_32, t_t, t_eqe, sigma=t_deqe)

plt.figure()
plt.errorbar(t_t, t_eqe, yerr=t_deqe, fmt='.', capsize=2)
plt.plot(t_t, t_32(t_t, *t_param,))
plt.xlabel('Temperature [Celsius]')
plt.ylabel('External quantum efficiency')
plt.grid(True)
plt.show()

#%%  TRESHOLD CURRENT

# Linear fit below and above the threshold
param_below, covm_below = opt.curve_fit(retta, current[mask_below], tot_p[mask_below], sigma= tot_p[mask_below]*0.05, absolute_sigma = True)
param_above, covm_above = opt.curve_fit(retta, current[mask], tot_p[mask], sigma= tot_p[mask]*0.05, absolute_sigma = True)


###     FIRST METHOD: LINEAR FIT
i_th[DATASET-1, 0] = -param_above[1]/param_above[0]
perr = np.sqrt(np.diag(covm_above))
di_th[DATASET-1, 0] = i_th[DATASET-1, 0]*np.sqrt(((perr[0]/param_above[0])**2)+((perr[1]/param_above[1])**2)-(2*(covm_above[0,1])/param_above[1]/param_above[0]))


###     SECOND METHOD: TWO-SEGMENT FIT
i_th[DATASET-1, 1] = - (param_above[1]-param_below[1])/(param_above[0]-param_below[0])
sigma_up = np.sqrt(covm_above[1,1] + covm_below[1,1])
sigma_down = np.sqrt(covm_above[0,0] + covm_below[0,0])
di_th[DATASET-1, 1] = i_th[DATASET-1, 1]*np.sqrt((sigma_up/(-param_above[1]+param_below[1]))**2 + (sigma_down/(param_above[0]-param_below[0]))**2)

#%% DERIVATIVE COMPUTATION OF THRESHOLD CURRENT
# Smooth data points, consider only points masked around threshold

CURRENT_TH = np.array([[6,12],
                        [8,12],
                        [9,13],
                        [10,15],
                        [13, 17]])

mask_th = np.array([current[i] > CURRENT_TH[DATASET-1, 0]   for i in range(NUM)])  &  np.array([current[i] < CURRENT_TH[DATASET-1, 1]   for i in range(NUM)])#take only the points above the threshold

#exclude specific data points
if DATASET == 1:
    mask_th[[16,17,22,25,29]] = False 
elif DATASET == 2:
    mask_th[17] = False
elif DATASET == 3:
    mask_th[[18,23]] = False
elif DATASET == 4:
    mask_th[[23,24,25]] = False
elif DATASET == 5:
    mask_th[[16,17]] = False

window_size=3
def mov_average(series, windowsize):
    i=0
    smooth_series = []
    while i<len(series)-windowsize+1:
        thiswindow = series[i:i+windowsize]
        windowaverage=sum(thiswindow)/windowsize
        smooth_series.append(windowaverage)
        i += 1
    
    return np.array(smooth_series)

smooth_i = mov_average(current[mask_th], window_size)
smooth_p = mov_average(tot_p[mask_th],window_size)

first_deriv_current = (smooth_i[1:] + smooth_i[:-1])/2
first_deriv = smooth_p[1:]-smooth_p[:-1]

smooth_1i = mov_average(first_deriv_current, window_size)
smooth_1p = mov_average(first_deriv, window_size)

second_deriv_current = (smooth_1i[1:] + smooth_1i[:-1])/2
second_deriv = smooth_1p[1:]-smooth_1p[:-1]

max_idx = np.argmax(second_deriv)
max_current = second_deriv_current[max_idx]
i_th[DATASET-1, 2]=max_current

print('Second derivative method ith = ', max_current)


plt.figure()
plt.plot(smooth_1i, smooth_1p, '.-', label='Diff1')
plt.plot(second_deriv_current, second_deriv, '.-', label='Diff2')
plt.grid(True)
plt.legend()

x1 = second_deriv_current[max_idx-1]
x2 = second_deriv_current[max_idx]
x3 = second_deriv_current[max_idx+1]
y1 = second_deriv[max_idx-1]
y2 = second_deriv[max_idx]
y3 = second_deriv[max_idx+1]
# Three point paraboloidal fit when there are very few data
denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2))
b = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3))
c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3)

print('ith optimized parabola is = ', -b/(2*a) )
i_th[DATASET-1, 3]=(second_deriv_current[max_idx]+second_deriv_current[max_idx-1])/2

#%%
xxx = np.linspace(9.5, 11.25, 50)
plt.plot(xxx, -(a*xxx**2+b*xxx+c ))
plt.show()

#%% FComparison methods plot

fig, ([ax1, ax3], [ax2, ax4]) = plt.subplots(2,2)
fig.suptitle('T = 296.001 K')

ax1.errorbar(current, tot_p, tot_p*0.05, fmt='.', capsize=2, label='Experimental data')
ax1.plot(current[mask], retta(current[mask], *param_above,), label='Linear fit')
ax1.vlines(x=i_th[DATASET-1, 0], ymin=-0.25, ymax=1, label='Threshold position')
ax1.set_xlabel('Current [mA]')
ax1.set_ylabel('Optical power [mW]')
ax1.title.set_text('Linear fit')
ax1.legend()
ax1.grid(True)


ax2.errorbar(current, tot_p, tot_p*0.05, fmt='.', capsize=2, label='Experimental data')
ax2.plot(current[mask], retta(current[mask], *param_above,), label='Linear fit')
ax2.plot(current[mask_below], retta(current[mask_below], *param_below), label='Linear fit')
ax2.vlines(x=i_th[DATASET-1, 1], ymin=-0.25, ymax=1, label='Threshold position')
ax2.set_xlabel('Current [mA]')
ax2.set_ylabel('Optical power [mW]')
ax2.legend()
ax2.title.set_text('Two segment fit')
ax2.grid(True)

ax3.plot(smooth_1i, smooth_1p, '.-', label='First derivative')
ax3.plot(second_deriv_current, second_deriv, '.-', label='Second derivative')
ax3.vlines(x=second_deriv_current[max_idx], ymin=-0.025, ymax=0.025, label='Maximum of second derivative')
ax3.set_xlabel('Current [mA]')
ax3.set_ylabel('Derivative [a.u.]')
ax3.legend()
ax3.title.set_text('Second derivative')
ax3.grid(True)

ax4.plot(smooth_1i, smooth_1p, '.-', label='First derivative')
ax4.plot(second_deriv_current, second_deriv, '.-', label='Second derivative')
ax4.plot(xxx, - 0.016895303655984*xxx**2 + 0.35672984578161*xxx - 1.8756456523420, '-', label='Paraboloidal fit')
ax4.vlines(x=-b/(2*a) , ymin=-0.025, ymax=0.025, label='Paraboloidal maximum')
ax4.set_xlabel('Current [mA]')
ax4.set_ylabel('Derivative [a.u.]')
ax4.legend()
ax4.title.set_text('Paraboloidal fit')
ax4.grid(True)

plt.show()

# ###     THIRD METHOD: FIRST DERIVATIVE
# f_deriv = np.gradient(tot_p)
# s_deriv = np.gradient(f_deriv)
# first_deriv = (tot_p[1:] - tot_p[:-1])/(current[1:] - current[:-1])
# first_deriv_current = (current[1:] + current[:-1])/2

# ###     FOURTH METHOD: SECOND DERIVATIVE
# second_deriv = (first_deriv[1:] - first_deriv[:-1])/(first_deriv_current[1:] - first_deriv_current[:-1])
# second_deriv_current = (first_deriv_current[1:] + first_deriv_current[:-1])/2
# i_th4 = second_deriv_current[np.argmax(second_deriv)]

# # UNCERTAINTY OF DERIVATIVE (APPROXIMATED AS FINITE DIFFERENCE)
# d_ptot = dtot_p[:,1]-dtot_p[:,0]
# sigm_up = np.sqrt(d_ptot[1:]**2+d_ptot[:-1]**2)

# i_th[DATASET-1]=i_th2
# di_th[DATASET-1]=d_ith2


# #%% POWER VS CURRENT PLOT
# fig, axs = plt.subplots(2,2)

# axs[0,0].fill_between(current, tot_p - dtot_p[:,0], tot_p + dtot_p[:,1], alpha=0.4)
# axs[0,0].plot(current, tot_p, '.-')
# axs[0,0].plot(np.arange(8,30), param_above[0]*np.arange(8,30) + param_above[1], '.-')
# axs[0,0].grid(True)

# axs[0,1].fill_between(current, tot_p - dtot_p[:,0], tot_p + dtot_p[:,1], alpha=0.4)
# axs[0,1].plot(current, tot_p, '.-')
# axs[0,1].plot(np.arange(8,30), param_above[0]*np.arange(8,30) + param_above[1], '.-')
# axs[0,1].plot(np.arange(0,15), param_below[0]*np.arange(0,15) + param_below[1], '.-')
# axs[0,1].grid(True)

# axs[1,0].plot(current, tot_p, '.-')
# axs[1,0].plot(first_deriv_current, first_deriv, '.-')
# axs[1,0].grid(True)

# axs[1,1].plot(current, tot_p, '.-')
# axs[1,1].plot(second_deriv_current, second_deriv, '.-')
# axs[1,1].grid(True)

# plt.xlabel('Current [mA]')
# plt.ylabel('Optical Power [mW]')
# plt.show()




#%% Characteristic temperature : threshold current as function of temperature.

T = np.empty(5)
for j in range(5):
        T[j] = np.loadtxt(name.format(j+1, 0), skiprows=2, max_rows=1, unpack=True)

ith = np.array([9.5167 , 10.5324, 11.7177, 13.0366,  15.5214])

dith = np.array([0.2757,
        0.2309,
        0.2733,
        0.2904,
        0.4245])

plt.figure()
plt.errorbar(T, np.log(ith), yerr=dith/ith,fmt='.r',label='Experimental data')
plt.xlabel('Temperature [°C]')
plt.ylabel('$ln(i_{th})$')
plt.grid(True)
plt.show()

param, covm = opt.curve_fit(retta, T, np.log(ith), sigma=dith/ith, absolute_sigma=True)
x = np.arange(15,50, 0.1)

plt.plot(x, param[0]*x+param[1], '-b',label='FIT')
plt.legend()

print('Characteristic temperature T0 =', 1/param[0], '+/-' , np.sqrt(covm[0,0]/param[0]**4), 'K')


