

#include <Metro.h>
#include <AltSoftSerial.h>    // Arduino build environment requires this
#include <wavTrigger.h>

wavTrigger wTrig;             // Our WAV Trigger object

Metro gLedMetro(500);         // LED blink interval timer
Metro gSeqMetro(6000);        // Sequencer state machine interval timer

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(13,OUTPUT);
  pinMode(12,OUTPUT);
  pinMode(11,OUTPUT);


  // If the Arduino is powering the WAV Trigger, we should wait for the WAV
  //  Trigger to finish reset before trying to send commands.
  delay(1000);

  // WAV Trigger startup at 57600
  wTrig.start();  
  delay(10);

  // Send a stop-all command and reset the sample-rate offset, in case we have
  //  reset while the WAV Trigger was already playing.
  wTrig.stopAllTracks();   
  wTrig.samplerateOffset(0); 

  // Sound volume control
  wTrig.masterGain(0);
  
}
int t=1;

void loop() {
  // put your main code here, to run repeatedly:

  int input = Serial.read();

  if (input == 49){
    wTrig.trackPlaySolo(t);
    //digitalWrite(11,HIGH);
    
    while(true){
      
      int sound = analogRead(A0);
      //Serial.println(sound);
      
      if (sound > 11){
        
        digitalWrite(12,HIGH);
        delay(3000);
        digitalWrite(11,LOW);
        digitalWrite(12,LOW);
        t++;
        
        for (int i = 0; i <20; i++){ //30s
    
          digitalWrite(11,HIGH);
          delay(1000);
          digitalWrite(11,LOW);
          digitalWrite(12,HIGH);
          delay(1000);
          digitalWrite(12,LOW);
          digitalWrite(13, HIGH);
          delay(1000);
          digitalWrite(13,LOW);

        }
        //*/
      break;
      }
      
    }
  }
}
