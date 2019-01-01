#=============== SETUP ===============

# Import libraries
import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt
import csv


##=============== FUNCTIONS ===============

# Retrieve data from CSV
def GetData(fileName):
    list_time = []
    list_rawRed = []
    list_rawIR = []
    
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile,delimiter=',')

        for row in reader:
            list_time.append(int(row[0]))
            list_rawRed.append(int(row[1]))
            list_rawIR.append(int(row[2]))
    
    time = np.asarray(list_time)
    rawRed = np.asarray(list_rawRed)
    rawIR = np.asarray(list_rawIR)
    
    return time, rawRed, rawIR


# Plot N signals for comparison
def PlotSignals(x, *y, offset=0):

    N = len(y)
    plotKey = N*100 + 11
    
    plt.figure()
    for i in range(0,N):
        plt.subplot(plotKey)
        plt.plot(x[offset:], y[i][offset:])
        plotKey = plotKey+1
    plt.show()
    
    return


# Get the Fast Fourier Transform of a signal
def GetFFT(signal, sampleFreq, offset = 0):
    
    N = signal.size
    f = np.linspace(0, sampleFreq, N)
    frequencies = f[offset:(N//2)]
    
    fft = np.fft.rfft(signal)
    magnitudes = np.abs(fft)[offset:(N//2)] / N
    
    return frequencies, magnitudes


# Bandpass filter
def Bandpass(x, f1, f2, fs, order = 6):
    fc1 = f1 / (0.5*fs)
    fc2 = f2 / (0.5*fs)
    
    b,a = signal.butter(order, [fc1,fc2], btype='bandpass')
    y = signal.lfilter(b, a, x)
    
    return y


# Low pass filter
def Lowpass(x, fc, order=6):
    
    b,a = signal.butter(order, fc, btype='low')
    y = signal.lfilter(b, a, x)
    
    return y


# High pass filter
def Highpass(x, fc, order=6):
    
    b,a = signal.butter(order, fc, btype='high')
    y = signal.lfilter(b, a, x)
    
    return y


# Perform frequency domain processing on epoch
def FrequencyAnalysis(signal,prev,f1,f2,fs):
    
    # Perform FFT
    freq, mags = GetFFT(signal,fs,offset=0)
    
    # Get bounds of examination - assume HR changes by at most 4.8bpm/s
    lower = max(prev-0.4,f1)
    lower = min(f2-0.8,lower)
    upper = lower + 0.8
    
    # Select region of data
    lower = int(lower/(freq[1] - freq[0]))
    upper = int(upper/(freq[1] - freq[0]))
    freq = freq[lower:upper]
    mags = mags[lower:upper]
    
    # Find peaks
    peaks = []
    strengths = []
    
    for a in np.arange(1,len(mags)-1):
        
        # Check if it's a peak
        if mags[a-1] < mags[a] and mags[a+1] < mags[a]:
            peaks.append(freq[a])
            strengths.append(mags[a])
            
# =============================================================================
#     # Print peaks
#     m = (np.asarray(strengths))
#     m = m.astype(int)
#     print('mags:',m)
#     
#     HRs = (60*np.asarray(peaks))
#     HRs = HRs.astype(int)
#     print('peaks:',HRs)
# =============================================================================
    
    # Return the maximum value
    peak = freq[np.argmax(mags)]
    #PlotSignals(freq,mags)
    
    
    return peak, peaks



#=============== MAIN ===============

# Set up key variables
fileName = '2018-10-06_Rest_Wrist-F_01.csv'

fs = 12.5
f1 = 0.5
f2 = 3.5

initEpoch = 15
epoch = 10
repeat = 5
HR = [1.20]


# Import data and set up parsing
time, rawRed, rawIR = GetData(fileName)

start = int(np.floor(initEpoch*fs))
size = np.shape(rawIR)[0]
indices = np.arange(start,size,int(repeat*fs))


# Parse through dataset and perform segment analysis
for i in indices[0:20]:
    
    # Retrieve current segment, perform BPF, retrieve epoch
    IR = rawIR[(i-start):i]
    BPF_IR = Bandpass(IR,f1,f2,fs,2)[int((initEpoch-epoch)*fs/2):]

    # Perform frequency analysis on epoch, retrieve most prominent peaks    
    rate, peaks = FrequencyAnalysis(BPF_IR,HR[len(HR)-1],f1,f2,fs)
    HR.append((rate + HR[len(HR)-1]) / 2 )
    
    # Retrieve time data for epoch
    epochStart = (i-start)+int((initEpoch-epoch)*fs/2)
    t = time[epochStart:i]
    
    #PlotSignals(t,BPF_IR)
    
    #print('HR:',int(60*HR[len(HR)-1]))
    #print(' ')



HR = (60*np.asarray(HR))
HR = HR.astype(int)
print(HR)
#xvals = np.arange(0,np.shape(HR)[0])
#PlotSignals(xvals,HR)

# =============================================================================
# # Working on a subset of data
# subset = 10                            # interval number
# epoch = 10                             # seconds of data [s]
# 
# lower = int(np.floor(subset*epoch*fs))
# upper = int(np.floor((subset+1)*epoch*fs))
# 
# t = time[lower:upper]
# R = rawRed[lower:upper]
# IR = rawIR[lower:upper]
# 
# 
# BPF_R = Bandpass(R,f1,f2,fs,2)
# BPF_IR = Bandpass(IR,f1,f2,fs,2)
# PlotSignals(t,R,IR,BPF_R,BPF_IR,offset=63)
# 
# 
# # FFT on the subset of data
# freq, BPF_R_FFT = GetFFT(BPF_R[63:],fs,offset=0)
# freq, BPF_IR_FFT = GetFFT(BPF_IR[63:],fs,offset=0)
# PlotSignals(freq,BPF_R_FFT,BPF_IR_FFT)
# 
# =============================================================================
