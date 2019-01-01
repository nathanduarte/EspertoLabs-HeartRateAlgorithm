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

fileName = '2018-09-25_Data_02.csv'
f1 = 0.5 / 6.25
f2 = 3.5 / 6.25

time, rawRed, rawIR = GetData(fileName)


# Red Vals
noDC_Red = Highpass(rawRed,f1,10)

band_Red1 = Bandpass(rawRed,f1,f2,8)
band_Red2 = Bandpass(np.power(rawRed,2),f1,f2,8)
band_Red3 = np.power(Bandpass(rawRed,f1,f2,8),2)
#PlotSignals(time,band_Red1,band_Red2,band_Red3,offset=200)

freqs_Red, band_Red1_FFT = GetFFT(band_Red1,12.5)
freqs_Red, band_Red2_FFT = GetFFT(band_Red2,12.5)
freqs_Red, band_Red3_FFT = GetFFT(band_Red3,12.5)

f3 = 10*(freqs_Red[1] - freqs_Red[0])
band_Red2_FFT_LP = np.power(Lowpass(band_Red2_FFT,f3,8),10)

PlotSignals(freqs_Red,band_Red1_FFT,band_Red2_FFT,band_Red2_FFT_LP,band_Red3_FFT)


# IR Vals
noDC_IR = Highpass(rawIR,f1,10)
band_IR = Bandpass(rawIR,f1,f2,8)

freqs, rawIR_FFT = GetFFT(rawIR,12.5)
freqs, band_IR_FFT = GetFFT(band_IR,12.5)