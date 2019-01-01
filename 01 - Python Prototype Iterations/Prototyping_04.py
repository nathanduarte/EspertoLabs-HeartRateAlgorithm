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
    

# Bandpass filter
def Bandpass(x, f1, f2, order = 6):
    
    b,a = signal.butter(order, [f1,f2], btype='bandpass')
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

fileName = '2018-10-08_Rest_Wrist-F_01_NMD.csv'
fs = 12.5
f1 = 0.5 / (0.5*fs)
f2 = 3.5 / (0.5*fs)

time, rawRed, rawIR = GetData(fileName)

# Try magnifying the signal
rawRed = np.power(rawRed,2)
rawIR = np.power(rawIR,2)

# BPF on the signals
BPFRed = Bandpass(rawRed,f1,f2,4)
BPFIR = Bandpass(rawIR,f1,f2,4)
PlotSignals(time,BPFRed,BPFIR,offset=100)


# Get FFT of the raw signals
freqs_Red, rawRed_FFT = GetFFT(rawRed,fs)
freqs_IR, rawIR_FFT = GetFFT(rawIR,fs)

f3 = 20*(freqs_Red[1] - freqs_Red[0])
rawRed_FFT_LP = Lowpass(rawRed_FFT,f3,1)
rawIR_FFT_LP = Lowpass(rawIR_FFT,f3,1)

off = int(0.4/(freqs_Red[1] - freqs_Red[0]))
PlotSignals(freqs_Red,rawRed_FFT,rawIR_FFT,rawRed_FFT_LP,rawIR_FFT_LP,
            offset=off)
















