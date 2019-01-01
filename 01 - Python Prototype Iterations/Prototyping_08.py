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
def FrequencyAnalysis(signal,f1,f2,fs):
    
    # Perform FFT
    freq, mags = GetFFT(signal,fs,offset=0)
    #PlotSignals(freq,mags)
    
    # Select region of data
    lower = int(f1/(freq[1] - freq[0]))
    upper = int(f2/(freq[1] - freq[0]))
    freq = freq[lower:upper]
    mags = mags[lower:upper]
    
    # Find the prominent peaks
    intervalSize = 5
    a = int((intervalSize+1)/2)
    
    peakFreqs = []
    peakMags = []

    for i in np.arange(a,len(mags)-a):
        x = freq[i-a:i+a]
        y = mags[i-a:i+a]
        
        mean = np.mean(y)
        sd = np.std(y)
        
        if np.argmax(y) == a:
            if (y[a] - mean)/sd > 1.2:
                peakFreqs.append(x[a])
                peakMags.append(y[a])
    
    peakFreqs = np.asarray(peakFreqs)
    peakMags = np.asarray(peakMags)
    
    #PlotPoints(freq,mags,peakFreqs,peakMags)
    
    return peakFreqs


# Perform time domain processing on epoch
def TimeAnalysis(signal,times,frequencies):
    
    # Determine intervals
    intervals = (12.5/frequencies).astype(int)
    for i in np.arange(0,np.shape(intervals)[0]):
        if intervals[i] % 2 == 0:
            intervals[i] = intervals[i] - 1

    # Find peaks
    halves = ((intervals-1)/2).astype(int)
    
    #halves = np.arange(1,1+2*np.shape(frequencies)[0],2)
    meanPeakDiff = []
    
    for i in np.arange(0,np.shape(halves)[0]):
        
        bot = halves[i]
        top = int(np.shape(signal)[0] - 1 - bot)
        
        tVals = []
        xVals = []
        
        for j in np.arange(bot,top):
            t = times[j-bot:j+bot]
            x = signal[j-bot:j+bot]
            
            if np.argmax(x) == bot:
                tVals.append(t[bot])
                xVals.append(x[bot])
                
        tVals = np.asarray(tVals)
        xVals = np.asarray(xVals)
        
        #PlotPoints(times,signal,tVals,xVals)
        
        rawFirstDiffs = np.diff(tVals) / 1000
        avg = np.mean(rawFirstDiffs)
        sd = np.std(rawFirstDiffs)
        deviations = np.abs(rawFirstDiffs - avg)/sd
    
        firstDiffs = np.extract(deviations < 1.5, rawFirstDiffs)
        meanPeakDiff.append(np.mean(firstDiffs))
        
    peakFreqs = 1 / np.asarray(meanPeakDiff)
    errors = np.abs(frequencies - peakFreqs)
    rate = frequencies[np.argmin(errors)]
    
    print(frequencies)
    print(peakFreqs)
    print(int(rate*60))

    #PlotSignals(times,signal)
    #PlotPoints(times,signal,peakTimes,peakMags)
    
    return rate


#=============== MAIN ===============

# Set up key variables
fileName = '2018-10-22_Excercise_Wrist-F_01_ND.csv'

fs = 12.5
f1 = 0.5
f2 = 3.5

initEpoch = 18
epoch = 10
repeat = 5
HR = []


# Import data and set up parsing
time, rawRed, rawIR = GetData(fileName)
PlotSignals(time,rawIR)

start = int(np.floor(initEpoch*fs))
size = np.shape(rawIR)[0]
indices = np.arange(start,size,int(repeat*fs))


# Parse through dataset and perform segment analysis
for i in indices[50:]:
    
    # Retrieve current segment, perform BPF, retrieve epoch
    IR = rawIR[(i-start):i]
    BPF_IR = Bandpass(IR,f1,f2,fs,2)[int((initEpoch-epoch)*fs/2):]

    # Perform frequency analysis on epoch, retrieve most prominent peaks
    pkFreqs = FrequencyAnalysis(BPF_IR,f1,f2,fs)
    
    # Retrieve time data for epoch
    epochStart = (i-start)+int((initEpoch-epoch)*fs/2)
    t = time[epochStart:i]
    
    # Perform time analysis on epoch, determine a rate
    rate = TimeAnalysis(BPF_IR,t,pkFreqs)
    
    # Determine overall rate and add it to the list
    if len(HR) > 0:
        HR.append((rate + HR[len(HR)-1]) / 2 )
    else:
        HR.append(rate)


HR = (60*np.asarray(HR))
HR = HR.astype(int)
print(HR)

xvals = np.arange(0,np.shape(HR)[0]*5,5)
PlotSignals(xvals,HR)
