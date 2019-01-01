# Import libraries
import serial
import threading
import numpy as np
import matplotlib.pyplot as plt



#=============== SETUP ===============

# Set up serial stream
arduino = serial.Serial('COM4',115200,timeout=0.1)

# Set up FIFO windows
windowCap = 20
windowSize = 0
windowCurrent = 0
windowFirst = 0

time = np.arange(windowCap)
rawRed = np.arange(windowCap)
rawIR = np.arange(windowCap)



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


# Remove oldest data from window and add new data to window to window queue
def UpdateWindow():
    global windowCap
    global windowSize
    global windowCurrent
    global windowFirst
    
    global time
    global rawRed
    global rawIR
    
    # Remove oldest data
    if(windowSize != 0):
        windowFirst = (windowFirst+1)%windowCap
        windowSize = windowSize-1
        
    # Collect new data
    incomingVals = GetNewVals()

    # Add new data to window
    time[windowCurrent] = incomingVals[0]
    rawRed[windowCurrent] = incomingVals[1]
    rawIR[windowCurrent] = incomingVals[2]
    
    windowCurrent = (windowCurrent+1)%windowCap
    windowSize = windowSize+1


# Plotting the graphs
def plotValues(time, sig1, sig2):
    global windowCurrent
    global windowFirst
    
    if(windowFirst <= windowCurrent):
        t = time[windowFirst:windowCurrent:1]
        s1 = sig1[windowFirst:windowCurrent:1]
        s2 = sig2[windowFirst:windowCurrent:1]
    else:
        
        t = time[windowFirst:] + time[:windowCurrent]
        s1 = sig1[windowFirst:] + sig1[:windowCurrent]
        s2 = sig2[windowFirst:] + sig2[:windowCurrent]
    
    plt.subplot(211)
    plt.plot(t, s1)
    
    plt.subplot(212)
    plt.plot(t, s2)
    
    plt.show()
    
    return
    

# Bandpass filter
    
    

#=============== MAIN ===============

i = 0
while(i < 500):
    UpdateWindow()
    
    if(windowCurrent % 2 == 0):
        plotValues(time, rawRed, rawIR)
        
    print(time[windowCurrent], end='')
    print(',', end='')
    print(rawRed[windowCurrent], end='')
    print(',', end='')
    print(rawIR[windowCurrent])
    
    i = i+1

arduino.close()