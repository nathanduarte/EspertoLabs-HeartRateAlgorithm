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
    
    return


# Plot points on a signal
def PlotPoints(x, y, xpts, ypts, offset=0):

    plt.figure()
    plt.plot(x,y)
    plt.plot(xpts,ypts,'o')
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
    PlotSignals(freq,mags)
    
    # Get bounds of examination
    maxRateOfChange = 4.8
    lower = max(prev - (5*maxRateOfChange/60),f1)
    lower = min(f2 - (2*5*maxRateOfChange/60),lower)
    upper = lower + (2*5*maxRateOfChange/60)
    
    # Select region of data
    lower = int(lower/(freq[1] - freq[0]))
    upper = int(upper/(freq[1] - freq[0]))
    freq = freq[lower:upper]
    mags = mags[lower:upper]
    
    # Find the prominent peaks
    peaks = []
    strengths = []

    for a in np.arange(2,len(mags)-2):
        if mags[a-1] < mags[a] and mags[a+1] < mags[a]:
# =============================================================================
#            # More thorough check - doesn't always return a peak, though 
#            if mags[a-2] < mags[a] and mags[a+2] < mags[a]:
#                 if mags[a] > np.mean(mags[a-2:a+3]):
#                     peaks.append(freq[a])
#                     strengths.append(mags[a])
# =============================================================================
                    
            peaks.append(freq[a])
            strengths.append(mags[a])
    
    peaks = np.asarray(peaks)
    strengths = np.asarray(strengths)
    
    # Find the maximum value
    peak = freq[np.argmax(mags)]
    
    #PlotPoints(freq,mags,peaks,strengths)
    
    return peak, peaks, strengths


# Perform time domain processing on epoch
def TimeAnalysis(signal,times,frequencies,magnitudes):
    
    # Determine interval for peak detection
    weight_sum = np.sum(frequencies * magnitudes) / np.sum(magnitudes)
    scale = 1
    interval = int(scale*(1/weight_sum)*12.5)
    
    if interval % 2 == 0:
        interval = interval - 1
    #print(interval)
    
    # Find peaks
    peakTimes = []
    peakMags = []
    bot = int((interval-1)/2)
    top = int(np.shape(signal)[0] - 1 - bot)
    
    for i in np.arange(bot,top):
        t = times[i-bot:i+bot]
        x = signal[i-bot:i+bot]
        if np.argmax(x) == bot:
            peakTimes.append(t[np.argmax(x)])
            peakMags.append(x[np.argmax(x)])
    
    # Convert peak times and magnitudes to np array
    peakTimes = np.asarray(peakTimes)
    peakMags = np.asarray(peakMags)
    
    # Find first differences and remove differences with an sd > 2  
    rawFirstDiffs = np.diff(peakTimes) / 1000
    
    avg = np.mean(rawFirstDiffs)
    sd = np.std(rawFirstDiffs)
    deviations = np.abs(rawFirstDiffs - avg)/sd
    
    firstDiffs = np.extract(deviations < 1.5, rawFirstDiffs)
    rate = 1 / np.mean(firstDiffs)
    
    #PlotSignals(times,signal)
    PlotPoints(times,signal,peakTimes,peakMags)
    
    return rate


#=============== MAIN ===============

# Set up key variables
fileName = '2018-10-18_Excercise_Wrist-F_03_ND.csv'

fs = 12.5
f1 = 0.5
f2 = 3.5

initEpoch = 15
epoch = 10
repeat = 5
HR = [1.2]


# Import data and set up parsing
time, rawRed, rawIR = GetData(fileName)

start = int(np.floor(initEpoch*fs))
size = np.shape(rawIR)[0]
indices = np.arange(start,size,int(repeat*fs))


# Parse through dataset and perform segment analysis
for i in indices:
    
    # Retrieve current segment, perform BPF, retrieve epoch
    IR = rawRed[(i-start):i]
    BPF_IR = Bandpass(IR,f1,f2,fs,2)[int((initEpoch-epoch)*fs/2):]

    # Perform frequency analysis on epoch, retrieve most prominent peaks    
    freqRate, freqPeaks, freqMags = FrequencyAnalysis(BPF_IR,HR[len(HR)-1],f1,f2,fs)
    
    # Retrieve time data for epoch
    epochStart = (i-start)+int((initEpoch-epoch)*fs/2)
    t = time[epochStart:i]
    
    # Perform time analysis on epoch, determine a rate
    timeRate = TimeAnalysis(BPF_IR,t,freqPeaks,freqMags)
    
    # Determine overall rate and add it to the list
    rate = (freqRate + timeRate)/2 
    HR.append((rate + HR[len(HR)-1]) / 2 )


HR = (60*np.asarray(HR))
HR = HR.astype(int)
print(HR)

xvals = np.arange(0,np.shape(HR)[0]*5,5)
PlotSignals(xvals,HR)