#%%
#==============================================================================
# Imports
#==============================================================================
from math import *
import cmath
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import sys
import os
from PyQt4.QtGui import *

#Import custom modules
from physics import *

sns.set(font_scale=1.5)
sns.set_style("ticks")
sns.set_palette(palette='deep')
sns.set_color_codes(palette='deep')
mpl.rcParams.update({'font.family': 'serif', 'font.serif':'DejaVu Serif'})

#%%
#==============================================================================
# User inputs
#==============================================================================

#GUI prompt
#a = QApplication(sys.argv) 
#w = QWidget()
#w.resize(320, 240)
#w.setWindowTitle("Select Directory")
#Directory = QFileDialog.getExistingDirectory(w, 'Select Directory', os.curdir)

#Or comment out the above and hard code the directory in:
Directory = os.curdir + '/Reflectivity_Data/2016-06-01/run6_-0.5to0.5psFine_lowestEnergy'

binsize = 0.1 #Size of bins (Volts)
minbin = 1.4
maxbin = 1.9

CalibrationSlope = 1.0
Calibration_yIntercept = 0

SelectPumpFluence = 1.8 #Desired pump fluence at which a lineout is taken (J/cm^2)
SelectDelay = -10 #Desired delay at which a lineout is taken (ps)

Connect_Points = False #Choose whether or not to connect the data points with lines

#%%
#==============================================================================
# Create list of file names for positive and negative delays
#==============================================================================
DelayPath = Directory + '/pos_*[!Offset]'
OffsetPath = Directory + '/pos_*Offset'

#Load in all of the filenames, sorted according to the digits at the end
DelayFilenames = sorted(glob.glob(DelayPath),key=lambda x: int(re.findall(r'_(\d*?)\Z',x)[0]))
OffsetFilenames = sorted(glob.glob(OffsetPath),key=lambda x: int(re.findall(r'_(\d*?)Offset',x)[0]))

#%%
#==============================================================================
# Read files and create arrays for signals
#==============================================================================
bins = np.arange(minbin,maxbin+binsize,binsize)

BinnedProbeVoltage = np.zeros((len(DelayFilenames),len(bins)))
MeanProbeVoltage = np.zeros(len(DelayFilenames))
MeanPumpVoltage = np.zeros(len(DelayFilenames))
MeanOfRatios = np.zeros(len(DelayFilenames))
Delay = np.zeros(len(DelayFilenames))

for i,filename in enumerate(DelayFilenames):
    Delay[i] = float(filename[(len(Directory)+5):(len(Directory)+13)])*-2/0.3 #ps of delay
        
    Data = pd.read_csv(filename, delimiter=",").values
    MeanPumpVoltage[i] = np.mean(Data[:,1],axis=0) #Input is inverted
    MeanProbeVoltage[i] = np.mean(Data[:,0],axis=0)
    MeanOfRatios[i] = np.sum(Data[:,0]/Data[:,1])/len(Data[:,1])
    
    #Energy binning:
    BinnedPumpVoltage = np.digitize(-Data[:,1],bins)
    for j in range(len(bins)):
        indices = np.where(BinnedPumpVoltage==j)
        BinnedProbeVoltage[i,j] = np.mean(Data[indices,0]/Data[indices,1])

RatioOfMeans = MeanProbeVoltage/MeanPumpVoltage

#Repeat for the offsets
OffsetBinnedProbeVoltage = np.zeros((len(OffsetFilenames),len(bins)))
OffsetMeanProbeVoltage = np.zeros(len(OffsetFilenames))
OffsetMeanPumpVoltage = np.zeros(len(OffsetFilenames))
OffsetMeanOfRatios = np.zeros(len(OffsetFilenames))
OffsetDelay = np.zeros(len(OffsetFilenames))

for i,filename in enumerate(OffsetFilenames):
    OffsetDelay[i] = float(filename[(len(Directory)+5):(len(Directory)+13)])*2/0.3 #ps of delay
    
    Data = pd.read_csv(filename, delimiter=",").values
    OffsetMeanPumpVoltage[i] = np.mean(Data[:,1],axis=0) #Input is inverted
    OffsetMeanProbeVoltage[i] = np.mean(Data[:,0],axis=0)
    OffsetMeanOfRatios[i] = np.sum(Data[:,0]/Data[:,1])/len(Data[:,1])
    
    #Energy binning:
    OffsetBinnedPumpVoltage = np.digitize(-Data[:,1],bins)
    for j in range(len(bins)):
        indices = np.where(OffsetBinnedPumpVoltage==j)[0]
        OffsetBinnedProbeVoltage[i,j] = np.mean(Data[indices,0]/Data[indices,1])
    
OffsetRatioOfMeans = OffsetMeanProbeVoltage/OffsetMeanPumpVoltage

#%%
#==============================================================================
# Subtract each delay by its associated offset
#==============================================================================
MeanProbeVoltage = MeanProbeVoltage-OffsetMeanProbeVoltage
MeanPumpVoltage = MeanPumpVoltage-OffsetMeanPumpVoltage
MeanOfRatios = MeanOfRatios-OffsetMeanOfRatios
RatioOfMeans = RatioOfMeans-OffsetRatioOfMeans
BinnedProbeVoltage = BinnedProbeVoltage - OffsetBinnedProbeVoltage

#%%
#==============================================================================
# Convert pump voltage into fluence
#==============================================================================
SpotSizeFWHM = sqrt(1/pi*200*log(2)) #FWHM focal spot size in um
w0 = SpotSizeFWHM/sqrt(2*log(2)) #Beam waist radius in um

MeanPumpEnergy = MeanPumpVoltage*CalibrationSlope + Calibration_yIntercept
PumpEnergy = bins*CalibrationSlope + Calibration_yIntercept #in uJ
PumpFluence = PumpEnergy/(pi*w0**2)*100 #Average (not peak) fluence in J/cm^2

#%%
#==============================================================================
# Select a particular bin and delay for creating lineouts
#==============================================================================

#Find index at which the array values are closest to the desired lineout value
SelectPumpFluenceIndex = min(range(len(PumpFluence)),key=lambda i: abs(PumpFluence[i]-SelectPumpFluence))
SelectDelayIndex = min(range(len(Delay)),key=lambda i: abs(Delay[i]-SelectDelay))

#Closest matching values
SelectPumpFluence = PumpFluence[SelectPumpFluenceIndex] 
SelectDelay = Delay[SelectDelayIndex]

#%%
#==============================================================================
# Create contour plot of probe signal with lineouts
#==============================================================================
if Connect_Points is True:
    linestyle = '-'
else:
    linestyle = 'None'

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(12,10))
fig.delaxes(ax4)

#Range of delays to plot
minDelay = np.min(Delay)
maxDelay = np.max(Delay)

#Create arrays necessary for mesh plot (2D arrays for energy and delay, masked array for probe signal)
BinnedProbeVoltageMasked = np.ma.masked_where(np.isnan(BinnedProbeVoltage),BinnedProbeVoltage)[:,1:] #Mask the NaN and ignore the first bin
EnergyGrid,DelayGrid = np.meshgrid(PumpEnergy[1:],Delay) #Ignore the first bin

#Contour plot of reflectance
mesh = ax1.pcolormesh(DelayGrid, EnergyGrid, BinnedProbeVoltageMasked,cmap='viridis')
ax1.axvline(SelectDelay,linestyle='--',color='k')
ax1.axhline(SelectPumpFluence,linestyle='--',color='r')
ax1.set_xlim(minDelay,maxDelay)
ax1.set_xlabel('Delay (ps)')
ax1.set_ylabel('Pump Fluence ($\mathrm{J/cm^2}$)')
#mesh.set_clim(0,0.1)
plt.colorbar(mesh,ax=ax1,label='$\Delta R$ (V)')

#Lineout across constant delay
ax2.plot(EnergyGrid[0,:],BinnedProbeVoltageMasked[SelectDelayIndex,:],linestyle=linestyle,marker='o',color='k')
ax2.set_xlabel('Pump Fluence ($\mathrm{J/cm^2}$)')

#Lineout across constant pump fluence
ax3.plot(Delay,BinnedProbeVoltage[:,SelectPumpFluenceIndex],linestyle=linestyle,marker='.',color='r')
ax3.set_xlim(np.min(Delay),np.max(Delay))
ax3.set_xlabel('Delay (ps)')
ax3.set_ylabel('$\Delta R$ (V)')

plt.tight_layout()
#plt.savefig('Reflectivity_Plot.png')

#%%
#==============================================================================
# Create plot of average ratios and pump and probe voltages
#==============================================================================

fig2 = plt.figure(figsize=(12,8))
plt.subplot(221)
plt.plot(Delay,MeanOfRatios,linestyle=linestyle,marker='.',label='Ratio')
plt.title('Mean of Ratios')
plt.ylabel('Voltage Ratio (Probe/Pump)')

plt.subplot(222)
plt.plot(Delay,RatioOfMeans,linestyle=linestyle,marker='.')
plt.title('Ratio of Means')

plt.subplot(223)
plt.plot(Delay,MeanProbeVoltage,linestyle=linestyle,marker='.')
plt.title('Probe Voltage')
plt.xlabel('Delay (ps)')
plt.ylabel('Voltage (V)')

plt.subplot(224)
plt.plot(Delay,MeanPumpVoltage,linestyle=linestyle,marker='.')
plt.title('Pump Voltage')
plt.xlabel('Delay (ps)')

#plt.xlim(-2.0,2.0)
plt.legend(frameon=True)
plt.tight_layout()
plt.show()