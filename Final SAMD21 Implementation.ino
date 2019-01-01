#include <Wire.h>
#include "MAX30105.h"
extern "C" char *sbrk(int i);


//========================= SETUP =========================
//--------------- MAX3010x SETUP ---------------
MAX30105 particleSensor;
#define debug SerialUSB      // SerialUSB if SAMD21


//--------------- COMPLEX MATH SETUP ---------------
struct cplx
{
  double real;
  double imag;
};
typedef struct cplx Complex;


//--------------- INITIALIZE VARIABLES ---------------
#define EPOCH 10          // time duration in consideration
#define RATE_SIZE 5       // number of heart rate values being stored
#define SAMPLE_RATE 12.5  // sample rate in Hz

// Heart rate tracking variables
double heartRate[RATE_SIZE];      // storing heart rate values in Hz
double heartRateAverage = 0;      // average of past 5 heart rates
int heartRateTracker = 0;         // for initial tracking of heart rate for averages

// Filtering and storage variables
const int cap = (int)(EPOCH * SAMPLE_RATE);
int current = 0;
int filterTracker = 0;        // for initial tracking of amount of samples
double MA_IR[cap];            // IR value after moving average filter
double IR[cap];               // IR value after bandpass filter
double IIRcoeff[] = {2.2898224619, -2.0425692391, 1.2453912738, -0.7801010845, 0.2707568350, -0.0105030013};
double gain = 6.60452365;

// Frequency analysis variables
double pi = 3.1415962;
Complex cplxIR[cap];          // for FFT
double magnitudes[cap/2];     // for FFT
double frequencies[cap/2];    // for FFT
double f1 = 0.5;
double f2 = 3.5;
const int maxNumberOfPeaks = 5;       // maximum peaks that can be drawn from frequency analysis
double peakFreq[maxNumberOfPeaks];    // for peaks from frequency analysis
double peakMags[maxNumberOfPeaks];    // for peaks from frequency analysis
int numberOfPeaks = 0;                // counting how many peaks were found



//========================= SETUP LOOP =========================
void setup()
{
  //--------------- MAX3010x SETUP ---------------
  debug.begin(115200);
  //debug.println("Initializing MAX3010x...");

  if (!particleSensor.begin(Wire, I2C_SPEED_FAST)) // use default I2C port, 400kHz speed
  {
    debug.println("MAX30105 was not found. Please check wiring/power. ");
    while (1);
  }
  particleSensor.setup(); // Configure sensor. Use 6.4mA for LED drive


  //--------------- INITIALIZE COLLECTION --------------
  // Collect 1 epoch of data to initialize bandpass filter
  for( int i = 0; i < cap; i++ )
  {
    AddNewValue(particleSensor.getIR());
  }
  
}



//========================= MAIN ALGORITHM =========================
void loop()
{
  //--------------- COLLECT DATA --------------
  AddNewValue(particleSensor.getIR());
  //debug.println(IR[current]);


  //--------------- PERFORM ANALYSIS --------------
  // Check for epoch
  if(current == 0){

    // UPDATE HR AVERAGE
    if(heartRateTracker >= RATE_SIZE){
      heartRateAverage = findMean(heartRate,RATE_SIZE);
      debug.print("HR: ");
      debug.println((int)(heartRateAverage*60));
    }


    // FREQUENCY ANALYSIS
    FrequencyAnalysis(IR);


    // TIME ANALYSIS
    double timeRate = 0;
    if(numberOfPeaks != 0){
      if(numberOfPeaks == 1){
        AddHeartRate(peakFreq[0]);
      }
      else{
        timeRate = TimeAnalysis(IR,cap);
        double differences[numberOfPeaks];
        int lowestDifferenceIndex = 0;
        for(int i = 0; i < numberOfPeaks; i++){
          differences[i] = abs(peakFreq[i] - timeRate);
          if(differences[i] < differences[lowestDifferenceIndex]){
            lowestDifferenceIndex = i;
          }
        }
        AddHeartRate(peakFreq[lowestDifferenceIndex]);
      }
    }
    debug.print("Free RAM: ");
    debug.println(FreeRAM());
  }
}



//========================= FUNCTIONS =========================

// Apply a moving average and a bandpass filter to a new value
void AddNewValue(int newVal)
{
  if(filterTracker < 2){
    MA_IR[current] = newVal;
    IR[current] = MA_IR[current] / gain;
    filterTracker++;
  }
  else if(filterTracker < 6){
    int a = current + cap;
    MA_IR[current] = 0.2*newVal + 0.5*MA_IR[(a-1)%cap] + 0.3*MA_IR[(a-2)%cap];
    IR[current] = MA_IR[current] / gain;
    filterTracker++;
  }
  else{
    int a = current + cap;
    MA_IR[current] = 0.2*newVal + 0.5*MA_IR[(a-1)%cap] + 0.3*MA_IR[(a-2)%cap];
    IR[current] = IIRcoeff[0]*IR[(a-1)%cap] + IIRcoeff[1]*IR[(a-2)%cap] + IIRcoeff[2]*IR[(a-3)%cap] + IIRcoeff[3]*IR[(a-4)%cap] + IIRcoeff[4]*IR[(a-5)%cap] + 
      IIRcoeff[5]*IR[(a-6)%cap] + ((MA_IR[current] - 3*MA_IR[(a-2)%cap] + 3*MA_IR[(a-4)%cap] - MA_IR[(a-6)%cap]) / gain);
  }
  
  current = (current + 1) % cap;
  return;
}


// Apply a Fast Fourier Transform using the FFT helper function
void FFT(double mags[], double sig[], int n)
{
  Complex complexSig[n];
  
  for( int i = 0; i < n; i++ ){
    complexSig[i].real = sig[i];
    complexSig[i].imag = 0;
  }
  
  FFTHelper(mags, complexSig, n, 1);

  double a;
  for( int i = 0; i < n; i++ ){
    a = mags[i] / n;
    mags[i] = a;
  }
  return;
}


// A Fast Fourier Transform helper function
void FFTHelper(double mags[], Complex complexSig[], int n, int jump)
{
  if(jump < n){
    FFTHelper(mags, complexSig, n, jump * 2);
    FFTHelper(mags + jump, complexSig + jump, n, jump * 2);

    for( int i = 0; i < n; i += 2 * jump ) {
      Complex t;
      t.real = cos(-1 * pi * i / n)*(complexSig[i + jump].real);
      t.imag = sin(-1 * pi * i / n)*(complexSig[i + jump].real);

      double real1 = complexSig[i].real + t.real;
      double real2 = complexSig[i].real - t.real;
      double imag1 = complexSig[i].imag + t.imag;
      double imag2 = complexSig[i].imag - t.imag;
      
      mags[i / 2] = sqrt(real1*real1 + imag1*imag1);
      mags[(i + n)/2] = sqrt(real2*real2 + imag2*imag2);
    }
  }
  return;
}


// Retrieve frequencies for FFT
void GetFrequencies(double freqs[], double fs, double N)
{
  double interval = fs / (N * 2);
  freqs[0] = 0;
  for( int i = 1; i < N; i++ ){
    freqs[i] = freqs[i-1] + interval;
  }
  return;
}


// Perform frequency analysis
void FrequencyAnalysis(double sig[])
{
  double HR_firstValue = heartRate[0];        // save value because it gets lost in FFT
  
  // Perform Fast Fourier Transform
  FFT(magnitudes,sig,cap/2);
  GetFrequencies(frequencies,SAMPLE_RATE,cap/2);

  heartRate[0] = HR_firstValue;               // replace saved value

  // Select region of data for frequency analysis
  int lower;
  int upper;
  int intervalSize;
  double sdThresh;
  
  GetRegionOfData(&lower,&upper,&intervalSize,&sdThresh);

  double freq[upper-lower];
  double mags[upper-lower];
  for(int i = 0; i < upper-lower; i++){
    freq[i] = frequencies[i+lower];
    mags[i] = magnitudes[i+lower];
  }

  // Reset peak holder variables
  numberOfPeaks = 0;
  for(int i = 0; i < maxNumberOfPeaks; i++){
    peakFreq[i] = 0;
    peakMags[i] = 0;
  }
  
  // Find the peaks  
  int a = (int)((intervalSize-1)/2);
  double segmentFreq[intervalSize];
  double segmentMags[intervalSize];

  for(int i = 0; i < (upper-lower+1) - intervalSize; i++){
    
    // Isolate segment
    for(int j = i; j < i + intervalSize; j++){
      segmentFreq[j-i] = freq[j];
      segmentMags[j-i] = mags[j];
    }

    // Get key values
    double mean = findMean(segmentMags,intervalSize);
    double stdev = findStDev(segmentMags,mean,intervalSize);
    int maxIndex = findMaxIndex(segmentMags,intervalSize);

    // Add peak if it is sufficient
    if(maxIndex == a){
      if((segmentMags[maxIndex] - mean)/stdev > sdThresh){

        if(numberOfPeaks < maxNumberOfPeaks){
          peakFreq[numberOfPeaks] = segmentFreq[maxIndex];
          peakMags[numberOfPeaks] = segmentMags[maxIndex];
          numberOfPeaks++;
        }
        else{
          int lowestIndex = 0;
          for(int k = 0; k < maxNumberOfPeaks; k++){
            if(peakMags[k] < peakMags[lowestIndex]){
              lowestIndex = k;
            }
          }
          peakFreq[lowestIndex] = segmentFreq[maxIndex];
          peakMags[lowestIndex] = segmentMags[maxIndex];
        }
      }
    }
  }

  // Condition where no peaks are found
  if(numberOfPeaks == 0){
    peakFreq[0] = heartRateAverage;
    peakMags[0] = 1;
    numberOfPeaks++;
  }
  
  return;
}


// Perform time analysis
double TimeAnalysis(double sig[], int n)
{
  // Retrieve weighted sum
  double weightSum = 0;
  double denominator = 0;
  for(int i = 0; i < numberOfPeaks; i++){
    weightSum += peakFreq[i]*peakMags[i];
    denominator += peakMags[i];
  }
  weightSum = weightSum / denominator;
  
  // Get interval
  int interval = (int)(SAMPLE_RATE/weightSum);
  if(interval % 2 == 0){
    interval = interval - 1;
  }
  int a = (int)((interval-1)/2);

  // Set up peak detection
  int firstDiffsCap = 100;
  double rawFirstDiffs[firstDiffsCap];
  double firstDiffs[firstDiffsCap];
  int numOfFirstDiffs = 0;
  int lastIndex = 0;

  // Perform peak detection
  double segmentSignal[interval];
  for(int i = 0; i < (n + 1) - interval; i++){
    
    // Isolate segment
    for(int j = i; j < i + interval; j++){
      segmentSignal[j-i] = sig[j];
    }

    // Check if peak is in the middle
    int maxIndex = findMaxIndex(segmentSignal,interval);
    if(maxIndex == a){
      if(lastIndex == 0){
        lastIndex = i + a;
      }
      else{
        if(numOfFirstDiffs < firstDiffsCap){
          int currentIndex = i + a;
          rawFirstDiffs[numOfFirstDiffs] = (currentIndex - lastIndex) / SAMPLE_RATE;
          numOfFirstDiffs++;
          lastIndex = currentIndex;
        }
      }
    }
  }

  // Remove difference outliers
  double mean = findMean(rawFirstDiffs,interval);
  double stdev = findStDev(rawFirstDiffs,mean,interval); 

  int firstDiffCounter = 0;
  for(int i = 0; i < numOfFirstDiffs; i++){
    double deviation = (rawFirstDiffs[i] - mean)/stdev;
    if(deviation < 1.5){
      firstDiffs[firstDiffCounter] = rawFirstDiffs[i];
      firstDiffCounter++;
    }
  }

  double rate = 1 / findMean(firstDiffs,firstDiffCounter);
  
  return rate;
}


// For selecting a segment of the frequency spectrum
void GetRegionOfData(int* lower, int* upper, int* intervalSize, double* sdThresh)
{
  if(heartRateTracker < RATE_SIZE){
    (*lower) = (int)(f1/(frequencies[1] - frequencies[0]));
    (*upper) = (int)(f2/(frequencies[1] - frequencies[0]));
    (*intervalSize) = 7;
    (*sdThresh) = 1.0;
  }
  else{
    // Maximum rate of change
    double maxFreqPerFiveSec = 1.2;
    double low = max(heartRateAverage-(maxFreqPerFiveSec/2),f1);
    low = min(f2-maxFreqPerFiveSec,low);
    double high = low + maxFreqPerFiveSec;
    
    (*lower) = int(low/(frequencies[1] - frequencies[0]));
    (*upper) = int(high/(frequencies[1] - frequencies[0]));
    (*intervalSize) = 5;
    (*sdThresh) = 0.8;
  }
  return;
}


// Mean of a set
double findMean(double set[], int n)
{
  double avg = 0;
  for(int i = 0; i < n; i++){
    avg = avg + set[i];
  }
  return avg / n;
}


// Standard deviation of a set
double findStDev(double set[], double mean, int n)
{
  double stdev = 0;
  for(int i = 0; i < n; i++){
    stdev += (set[i] - mean)*(set[i] - mean);
  }
  return sqrt(stdev / n);
}


// Index of the maximum value of a set
int findMaxIndex(double set[], int n)
{
  int maxIndex = 0;
  for(int i = 0; i < n; i++){
    if(set[i] > set[maxIndex]){
      maxIndex = i;
    }
  }
  return maxIndex;
}


void AddHeartRate(double rate)
{
  if(heartRateTracker < RATE_SIZE){
    heartRate[heartRateTracker] = rate;
    debug.print("HR: ");
    debug.println((int)(60*rate));
    heartRateTracker++;
  }
  else{
    for(int i = 0; i < RATE_SIZE-1; i++){
      heartRate[i] = heartRate[i+1];
    }
    heartRate[RATE_SIZE - 1] = rate;
  }
  return;
}

int FreeRAM()
{
  char stack_dummy = 0;
  return &stack_dummy - sbrk(0);
}
