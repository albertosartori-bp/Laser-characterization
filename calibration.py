import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


#%%     DIFFERENT AVERAGES
temp = np.loadtxt("Calibration/Medie_{}.txt".format(0), skiprows=4, unpack=True)
average_current = 25 # mA
average_temp = 25
aver_wlength = np.empty((4, temp.shape[1]))
aver_dbm_power = np.empty((4, temp.shape[1]))

# read data
for i in range(4):
    aver_wlength[i,:], aver_dbm_power[i,:] = np.loadtxt("CalibrationA/Medie_{}.txt".format(i), skiprows=4, unpack=True)
aver_power = 10**(aver_dbm_power/10)

# plot data
plt.figure()
for i in range(4):
    plt.plot(aver_wlength[i], aver_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()

# plot differences with max-mean data
plt.figure()
for i in range(4):
    plt.plot(aver_wlength[i], aver_dbm_power[3] - aver_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()

# total power
aver_tot_p = np.zeros(4)
for i in range(4):
    aver_tot_p[i] = np.sum(aver_power[i])
print(aver_tot_p)





#     DIFFERENT RESOLUTION
#%% FIRST CASE: T = 45, I = 29 mA, Medium mean
temp = np.loadtxt("Calibration/Res_{}.txt".format(0), skiprows=4, unpack=True)
#res1_resolutions = np.array([0.5,1,5]) # nm
res1_wlength = np.empty((3, temp.shape[1]))
res1_dbm_power = np.empty((3, temp.shape[1]))

# read data
for i in range(3):
    res1_wlength[i,:], res1_dbm_power[i,:] = np.loadtxt("Calibration/Res_{}.txt".format(i), skiprows=4, unpack=True)
res1_power = 10**(res1_dbm_power/10)

# plot data
plt.figure()
for i in range(3):
    plt.plot(res1_wlength[i], res1_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()

# total power
res1_tot_p = np.zeros(3)
res1_tot_p = np.sum(res1_power, axis=1)
print(res1_tot_p)


#%% SECOND CASE: T = 25, I = 20 mA, No mean
temp = np.loadtxt("Calibration/Reso_{}.txt".format(0), skiprows=4, unpack=True)
res2_resolutions = np.array([ 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 5]) # nm (the first one is Full resolution, not shown here)
res2_wlength = np.empty((9, temp.shape[1]))
res2_dbm_power = np.empty((9, temp.shape[1]))

# read data
for i in range(9):
    res2_wlength[i,:], res2_dbm_power[i,:] = np.loadtxt("Calibration/Reso_{}.txt".format(i), skiprows=4, unpack=True)
res2_power = 10**(res2_dbm_power/10)

# plot data
plt.figure()
for i in range(7):
    plt.plot(res2_wlength[i], res2_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()


# total power
res2_tot_p = np.zeros(9)
res2_tot_p_skip = np.zeros(8)
for i in range(8):
    j = 0
    while j < res2_power[i+1].size:
        res2_tot_p_skip[i] += res2_power[i+1, j]
        j += 2*int(res2_resolutions[i]/res2_resolutions[0])
res2_tot_p = np.sum(res2_power, axis=1)
print(res2_tot_p)
print(res2_tot_p_skip)







#   RANGE
#%% FIRST CASE: T = 25, I = 25 mA, No averages, Full resolution
temp = np.loadtxt("Calibration/Range0_{}.txt".format(0), skiprows=4, max_rows=20000, unpack=True)
# i range sono: 1450-1650, 1449-1649, 1448-1648, 1447-1647
range1_wlength = np.empty((4, temp.shape[1]))
range1_dbm_power = np.empty((4, temp.shape[1]))

# read data
for i in range(4):
    range1_wlength[i,:], range1_dbm_power[i,:] = np.loadtxt("Calibration/Range0_{}.txt".format(i), skiprows=4, max_rows=20000, unpack=True)
range1_power = 10**(range1_dbm_power/10)

# plot data
plt.figure()
for i in range(4):
    plt.plot(range1_wlength[i], range1_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()


# total power
range1_tot_p = np.zeros(4)
range1_tot_p = np.sum(range1_power, axis=1)
print(range1_tot_p)


#%% SECOND CASE: T = 25, I = 25 mA, No averages, 2 nm resolution
temp = np.loadtxt("Calibration/Range1_{}.txt".format(0), skiprows=4, max_rows=20000, unpack=True)
# i range sono: 1450-1650, 1449-1649, 1448-1648, 1447-1647
range2_wlength = np.empty((4, temp.shape[1]))
range2_dbm_power = np.empty((4, temp.shape[1]))

# read data
for i in range(4):
    range2_wlength[i,:], range2_dbm_power[i,:] = np.loadtxt("Calibration/Range1_{}.txt".format(i), skiprows=4, max_rows=20000, unpack=True)
range2_power = 10**(range2_dbm_power/10)

# plot data
plt.figure()
for i in range(4):
    plt.plot(range2_wlength[i], range2_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()


# total power
range2_tot_p = np.zeros(4)
range2_tot_p = np.sum(range1_power, axis=1)
print(range2_tot_p)





#   STABILITY
#%% FIRST CASE: T = 25, I = 20 mA, No averages, Full resolution
temp = np.loadtxt("Calibration/Stability_0_{}.txt".format(1), skiprows=4, unpack=True)
# i range sono: 1450-1650, 1449-1649, 1448-1648, 1447-1647
stab1_wlength = np.empty((6, temp.shape[1]))
stab1_dbm_power = np.empty((6, temp.shape[1]))

# read data
for i in range(6):
    stab1_wlength[i,:], stab1_dbm_power[i,:] = np.loadtxt("Calibration/Stability_0_{}.txt".format(i+1), skiprows=4, unpack=True)
stab1_power = 10**(stab1_dbm_power/10)

# plot data
plt.figure()
for i in range(6):
    plt.plot(stab1_wlength[i], stab1_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()


# total power
stab1_tot_p = np.zeros(6)
stab1_tot_p = np.sum(stab1_power, axis=1)
print(stab1_tot_p)

#%% SECOND CASE: T = 25, I = 25 mA, No averages, Full resolution
temp = np.loadtxt("Calibration/Stability_1_{}.txt".format(1), skiprows=4, unpack=True)
# i range sono: 1450-1650, 1449-1649, 1448-1648, 1447-1647
stab2_wlength = np.empty((5, temp.shape[1]))
stab2_dbm_power = np.empty((5, temp.shape[1]))

# read data
for i in range(5):
    stab2_wlength[i,:], stab2_dbm_power[i,:] = np.loadtxt("Calibration/Stability_1_{}.txt".format(i+1), skiprows=4, unpack=True)
stab2_power = 10**(stab2_dbm_power/10)

# plot data
plt.figure()
for i in range(5):
    plt.plot(stab2_wlength[i], stab2_dbm_power[i], label='{}'.format(i+1))
    #plt.plot(wlength[i, max_back[i]], dbm_power[i, max_back[i]], 'k.')

plt.legend(loc='upper right')
plt.xlabel('wavelength [nm]')
plt.ylabel('power [dBm]')
plt.grid(True)
plt.show()


# total power
stab2_tot_p = np.zeros(5)
stab2_tot_p = np.sum(stab2_power, axis=1)
print(stab2_tot_p)