# EspertoLabsHR_F18
This algorithm was developed for the Esperto Watch in Fall 2018. Esperto Labs (www.espertolabs.ca) is a student design team at the University of Waterloo.


**System Overview:** The Watch uses a MAX30102 pulse plethysmography (PPG) sensor which has an infrared (IR) and red channel. For this heart rate (HR) detection algorithm, only the IR channel was used. A SAMD21 microcontroller was used for real-time processing and was connected to the MAX30102 sensor via an I2C connection. The embedded code was implemented in an Arduino C environment.


**Algorithm Design:** The algorithm uses a combined frequency-domain and time-domain approach to calculate heart rate. Heart rate is determined every 10 seconds using the past 10 seconds of data with a sampling rate of 12.5 Hz.

*Filtering*: First, a third-order Butterworth bandpass filter is used to isolate frequencies between 0.5Hz and 3.5Hz.

*Frequency-Domain*: Next, a Fast Fourier Transform (Cooley-Tukey implementation) is used to identify prominent frequencies in the signal with the assumption that one of the frequencies corresponds to the heart rate. The slope-sign-change method coupled with some statistical thresholding is used to identify prominent frequencies.

*Time-Domain*: Then, a similar peak-finding algorithm is applied to the filtered IR data in the time-domain; the window used for this time-domain peak finding is determined by the prominent frequencies identified in the frequency-domain. The first differences between the peaks found in the time-domain are calculated, those with too great a deviation are rejected, and the remaining first differences are used to determine the time-domain HR.

*Final HR*: Finally, the time-domain HR is compared with the prominent frequencies identified in the frequency-domain and the frequency that is closest to the time-domain HR is considered the current HR.


**Next Steps:** There are a number of key parameters which greatly affect the accuracy of this algorithm especially ones related to peak-finding (in the time and frequency domain) and ones related to HR storage and reporting. These parameters must be optimized and further user testing must be carried out prior to full integration with the Watch.


**Repository Contents:** This repository contains all iterations of algorithm prototyping (in Python) and algorithm implementation (in C). Datasets used for algorithm development are not included here.
