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


# Circular FIFO queue class for stimulating windows
class CircularQueue:
    def __init__(self, capacity):
        self.queue = np.arange[capacity]
        self.cap = capacity
        self.size = 0
        self.first = 0
        self.current = 0
        return
    
    def Enqueue(val):
        if(self.size == self.cap):
            return
        else:
            self.queue[self.current] = val
            self.size = self.size+1
            self.current = (self.current+1)%self.cap
            return
        
    def Dequeue():
        if(self.size == 0):
            return 0
        else:
            temp = self.queue[self.first]
            self.size = self.size-1
            self.first = (self.first+1)%self.cap
            return temp


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


# FFT on parts of the signal
def FFTonEpoch(x,epoch,fs):
    samples = np.shape(x)[0]
    intervals = int(np.floor((samples/fs)/epoch))
    
    for i in np.arange(0,intervals):
        lower = int(np.floor(i*epoch*fs))
        upper = int(np.floor((i+1)*epoch*fs))
        freq, mag = GetFFT(x[lower:upper],fs,offset=2)
        PlotSignals(freq,mag)
    
    return


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



#=============== MAIN ===============

fileName = '2018-10-03_Rest_Wrist_01.csv'
fs = 12.5
f1 = 0.5
f2 = 3.5


time, rawRed, rawIR = GetData(fileName)
#PlotSignals(time,rawRed,rawIR,offset=0)


# Information for subset analysis
initEpoch = 15                          # seconds of data to filter [s]
epoch = 10                              # seconds of data to analyze [s]
low = int(np.floor(initEpoch*fs))
size = np.shape(rawIR)[0]


# Perform frequency domain processing on epoch
def FrequencyHR(signal,prev,f1,f2,fs):
    
    # Perform FFT
    freq, mags = GetFFT(signal,fs,offset=0)
    
    # Get bounds of examination - need to make this more robust
    lower = int(np.maximum(f1,prev-0.5)/(freq[1] - freq[0]))
    upper = int(np.minimum(f2,prev+0.5)/(freq[1] - freq[0]))
    freq = freq[lower:upper]
    mags = mags[lower:upper]
    
    # Find peaks
    sd = np.std(mags)
    mean = np.mean(mags)
    sdevs = (mags - mean)/sd
    peak = freq[np.argmax(sdevs)]
    print(int(60*peak))
    
    return peak


# Subset analysis
HR = [1.1]          # intial HR value

for i in np.arange(low,size,int((epoch/2)*fs)):
    IR = rawIR[(i-low):i]
    BPF_IR = Bandpass(IR,f1,f2,fs,2)[int(epoch*fs/2):]
    
    # Perform frequency domain analysis on BPF_IR
    rate = FrequencyHR(BPF_IR,HR[len(HR)-1],f1,f2,fs)
    HR.append(rate)

    # need a method of updating/modifying previous values


xval = np.arange(0,len(HR))
HR = Lowpass(HR,0.3,2)
PlotSignals(xval,HR)


# Working on a subset of data
subset = 3                            # interval number
initEpoch = 15                        # seconds of data [s]

lower = int(np.floor(subset*initEpoch*fs))
upper = int(np.floor((subset+1)*initEpoch*fs))

t = time[lower:upper]
R = rawRed[lower:upper]
IR = rawIR[lower:upper]


BPF_R = Bandpass(R,f1,f2,fs,2)
BPF_IR = Bandpass(IR,f1,f2,fs,2)
PlotSignals(t,R,IR,BPF_R,BPF_IR,offset=63)


# FFT on the subset of data
freq, BPF_R_FFT = GetFFT(BPF_R[63:],fs,offset=0)
freq, BPF_IR_FFT = GetFFT(BPF_IR[63:],fs,offset=0)
PlotSignals(freq,BPF_R_FFT,BPF_IR_FFT)