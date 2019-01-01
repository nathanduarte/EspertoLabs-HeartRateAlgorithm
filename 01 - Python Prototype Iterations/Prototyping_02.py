# Import libraries
import serial
import numpy as np
import matplotlib.pyplot as plt



#=============== SETUP ===============

# Set up serial stream
arduino = serial.Serial('COM4',115200,timeout=0.1)

# Set up data collection
capacity = 200
time = np.arange(capacity)
rawRed = np.arange(capacity)
rawIR = np.arange(capacity)



##=============== FUNCTIONS ===============

# Retrieve new values from the serial plotter
def GetNewVals():
    global arduino
    
    rawVals = np.arange(3)
    
    # Read and convert to string
    serialFeed = arduino.readline()
    commaSeparated = serialFeed[0:len(serialFeed)-2]    
    commaSeparatedStr = str(commaSeparated)
    
    # Process data and store it in rawVals
    if(len(commaSeparated) != 0):
        a = commaSeparatedStr.split(',')
        rawVals[0] = a[0][2:]
        rawVals[1] = a[1]
        rawVals[2] = a[2][:-1]
    
    return rawVals


# Collect all new data
def CollectData():
    global capacity, time, rawRed, rawIR
    
    for i in range(1,10):
        bufferVals = GetNewVals()
    
    for i  in np.arange(0,capacity):
        incomingVals = GetNewVals()
        
        time[i] = incomingVals[0]
        rawRed[i] = incomingVals[1]
        rawIR[i] = incomingVals[2]
        
    return


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
        


# Plotting the graphs
def PlotValues(time, sig1, sig2):
    
    plt.subplot(211)
    plt.plot(time[50:], sig1[50:])
    
    plt.subplot(212)
    plt.plot(time[50:], sig2[50:])

# =============================================================================
#     plt.subplot(211)
#     plt.plot(time, sig1)
#     
#     plt.subplot(212)
#     plt.plot(time, sig2)
# =============================================================================

    plt.show()
    
    return
    

# Bandpass filters assuming 12.5Hz sample rate; 0.5Hz to 3.5Hz
def BandpassFilter1(x):
    y = np.arange(x.size)
    G = 3.632746572
    xCoeffs = [1,0,-2,0,1]
    xCoeffs = [m * G for m in xCoeffs]
    yCoeffs = [-0.1725312505,0.3077544838,-0.7469155717,1.5242126723]

    for i in range(0,len(yCoeffs)):
        y[i] = x[i]
    
    for n in range(len(yCoeffs),y.size):
        sumX = 0
        sumY = 0
    
        for j in range(0,len(xCoeffs)):
            sumX = sumX + (xCoeffs[j] * x[n-j])
        for k in range(0,len(yCoeffs)):
            sumY = sumY + (yCoeffs[k] * y[n+k-len(yCoeffs)])
        
        y[n] = sumX + sumY
        
    return y

# Seems worse
def BandpassFilter2(x):
    G = 12.09580694
    y = np.arange(x.size)
    xCoeffs = [1,0,-4,0,6,0,-4,0,1]
    xCoeffs = [m * G for m in xCoeffs]
    yCoeffs = [-0.0181370770,0.0775148538,-0.4118091539,
               1.3206433598,-2.3639018292,3.2502425045,
               -3.9179486593,3.0550134982]

    for i in range(0,len(yCoeffs)):
        y[i] = x[i]
    
    for n in range(len(yCoeffs),y.size):
        sumX = 0
        sumY = 0
    
        for j in range(0,len(xCoeffs)):
            sumX = sumX + (xCoeffs[j] * x[n-j])
        for k in range(0,len(yCoeffs)):
            sumY = sumY + (yCoeffs[k] * y[n+k-len(yCoeffs)])
        
        y[n] = sumX + sumY
        
    return y


#=============== MAIN ===============

CollectData()

PlotValues(time, rawRed, rawIR)

PlotValues(time, BandpassFilter1(rawRed), BandpassFilter1(rawIR))

#PlotValues(time, BandpassFilter2(rawRed), BandpassFilter2(rawIR))
    
arduino.close()