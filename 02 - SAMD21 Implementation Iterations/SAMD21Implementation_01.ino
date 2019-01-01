#include <Wire.h>
#include "MAX30105.h"



//========================= SETUP =========================
//--------------- MAX3010x SETUP ---------------
MAX30105 particleSensor;
#define debug SerialUSB      // SerialUSB if SAMD21


//--------------- INITIALIZE VARIABLES ---------------
#define EPOCH 10          // time duration in consideration
#define RATE_SIZE 10      // number of heart rate values being stored
#define SAMPLE_RATE 12.5  // sample rate in Hz

// Heart rate tracking variables
double heartRate[RATE_SIZE];      // storing heart rate values in Hz
double lastRatesAvg;              // to be fed into frequency analysis

// Filtering and storage variables
const int cap = (int)(EPOCH * SAMPLE_RATE);
int current = 0;
int tracker = 0;             // for initial tracking of amount of samples
double MA_IR[cap];           // IR value after moving average filter
double IR[cap];              // IR value after bandpass filter
double IIRcoeff[] = {2.2898224619, -2.0425692391, 1.2453912738, -0.7801010845, 0.2707568350, -0.0105030013};
double gain = 6.60452365;



//========================= SETUP LOOP =========================
void setup()
{
  //--------------- MAX3010x SETUP ---------------
  debug.begin(115200);
  debug.println("Initializing MAX3010x...");

  if (!particleSensor.begin(Wire, I2C_SPEED_FAST)) // Use default I2C port, 400kHz speed
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
    //debug.print(millis());
    //debug.print(",");
    debug.println(IR[i]);
  }
}



//========================= MAIN ALGORITHM =========================
void loop()
{
  AddNewValue(particleSensor.getIR());
  //debug.print(millis());
  //debug.print(",");
  debug.println(IR[current]);
}



//========================= FUNCTIONS =========================

// Apply a moving average and a bandpass filter to a new value
void AddNewValue(int newVal)
{
  if(tracker < 2){
    MA_IR[current] = newVal;
    IR[current] = MA_IR[current] / gain;
    tracker++;
  }
  else if(tracker < 6){
    int a = current + cap;
    MA_IR[current] = 0.2*newVal + 0.5*MA_IR[(a-1)%cap] + 0.3*MA_IR[(a-2)%cap];
    IR[current] = MA_IR[current] / gain;
    tracker++;
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
