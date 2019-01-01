#include "MAX30105.h"
//#include "complex.h"



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
#define RATE_SIZE 10      // number of heart rate values being stored
#define SAMPLE_RATE 12.5  // sample rate in Hz

// Heart rate tracking variables
double heartRate[RATE_SIZE];      // storing heart rate values in Hz
double lastRatesAvg;              // to be fed into frequency analysis

// Filtering and storage variables
const int cap = (int)(EPOCH * SAMPLE_RATE) + 3;   // add 3 here so we get a power of 2 for the FFT
int current = 0;
int tracker = 0;             // for initial tracking of amount of samples
double MA_IR[cap];           // IR value after moving average filter
double IR[cap];              // IR value after bandpass filter
double IIRcoeff[] = {2.2898224619, -2.0425692391, 1.2453912738, -0.7801010845, 0.2707568350, -0.0105030013};
double gain = 6.60452365;

// Frequency analysis variables
double pi = 3.1415962;
Complex cplxIR[cap];
double f1 = 0.5;
double f2 = 3.5;
double magnitudes[cap/2];
double frequencies[cap/2];

// Test data artificial data
const int testN = 128;
double testData[] = {16412.23588,16413.44717,16414.02253,16415.94546,16416.92812,16417.20898,15966.31271,14935.45433,13500.22102,11758.83732,9769.26823,7666.764215,5594.40554,3641.647035,1877.394786,367.3509592,-849.1279959,-1763.071358,-2382.773846,-2727.148332,-2828.358228,-2733.539663,-2492.444214,-2149.764005,-1746.870733,-1319.163408,-888.4115589,-468.5079897,-89.40822428,209.1065507,415.0565532,549.2447117,627.4009168,652.8817611,642.1375234,608.5412765,544.5313286,447.1839482,333.2850267,224.799731,133.3170894,58.12464837,-3.179734717,-49.40910792,-83.30056532,-110.5967399,-128.9314052,-133.0000349,-127.8247387,-118.3344106,-101.8001689,-81.19802112,-63.49781294,-48.37007522,-32.51175615,-16.08590519,-0.890966182,10.8156078,18.04724822,22.02310126,22.86610595,22.26366152,22.43178529,20.23424243,14.90643786,10.35909245,9.498357367,15.18652387,24.42658034,29.03571205,24.51312444,10.80195563,-6.497008311,-22.50472927,-35.35302798,-38.32362989,-30.87796391,-22.87061621,-17.1859119,-10.41353612,-3.384098612,2.128145183,4.831495976,2.67236589,-0.448457543,2.959287055,10.55107133,13.63235148,11.10258216,7.776886445,5.216353895,2.987106197,2.767018419,4.152764158,3.920024789,1.067578607,-3.174318897,-6.06928294,-6.345882603,-7.471293664,-9.46662178,-7.383847768,-1.176722894,5.449478604,9.890120595,10.43350168,7.04510219,2.24921108,0.197022587,0.782801018,0.115392849,-1.074212745,-0.532365083,1.689945278,5.556520922,8.361929625,6.988779689,2.87446746,-0.670495693,-0.926261444,0.617894127,-1.8431946,-7.308639177,-8.988666962,-6.501632888,-3.298506188,-1.714037382,-4.493335777};



//========================= SETUP LOOP =========================
void setup()
{
  //--------------- MAX3010x SETUP ---------------
  debug.begin(115200);
  debug.println("Initializing MAX3010x...");  
}



//========================= MAIN ALGORITHM =========================
void loop()
{
  FFT(magnitudes,testData,testN/2);
  GetFrequencies(frequencies,SAMPLE_RATE,testN/2);
  
  debug.println("1");
  for(int i = 0; i < testN/2; i++){
//    //debug.print("f: ");
    debug.print(frequencies[i]);
    debug.print(",");
//    //debug.print(" m: ");
    debug.println(magnitudes[i]);
  }
}



//========================= FUNCTIONS =========================

// Apply a Fast Fourier Transform using the FFT helper function
void FFT(double mags[], double sig[], int n)
{
  Complex complexSig[n];
  
  for( int i = 0; i < n; i++ ){
    complexSig[i].real = sig[i];
    complexSig[i].imag = 0;
  }
  
  FFTHelper(mags, complexSig, n, 1);

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

//    for( int i = 0; i < n; i++ ){
//      double a = mags[i];
//      mags[i] = abs(a);
//    }
  }
}


// Retrieve frequencies for FFT
void GetFrequencies(double freqs[], double fs, double N)
{
  double interval = fs / (N * 2);
  freqs[0] = 0;
  for( int i = 1; i < N; i++ ){
    freqs[i] = freqs[i-1] + interval;
  }
}


/*
// It's unclear how useful these are
void _fft(Complex buf[], int n)
{
  Complex out[n];
  for (int i = 0; i < n; i++){
    out[i] = buf[i];
  }
  
  fft_helper(buf, out, n, 1);
}

void _fft_helper(Complex buf[], Complex out[], int n, int jump)
{
  if(jump < n) {
    fft_helper(out, buf, n, jump * 2);
    fft_helper(out + jump, buf + jump, n, jump * 2);
    
    for( int i = 0; i < n; i += 2 * jump ) {
      Complex holder(0,-1 * pi * i / n);
      Complex t = holder * out[i + jump];
      buf[i / 2] = out[i] + t;
      buf[(i + n)/2] = out[i] - t;
    }
  }
}

PI = atan2(1, 1) * 4;
cplx buf[] = {1, 1, 1, 1, 0, 0, 0, 0};

fft(buf, 8);


N = signal.size
f = np.linspace(0, sampleFreq, N)
frequencies = f[offset:(N//2)]
    
    fft = np.fft.rfft(signal)
    magnitudes = np.abs(fft)[offset:(N//2)] / N
 */
