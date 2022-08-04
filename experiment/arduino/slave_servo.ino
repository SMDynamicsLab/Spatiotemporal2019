// Código para el arduino esclavo
//Recibe la señal y el tamaño del escalón del servo

#include <Wire.h>
#include <Servo.h>

Servo myservo;  // create servo object to control a servo

void setup() {
  Wire.begin(8);                // join i2c bus with address #8
  Wire.onReceive(receiveEvent); // register event
  myservo.attach(9);            // attaches the servo on pin 9 
  myservo.write(57);
  Serial.begin(9600);           // start serial for output
}

void loop() {
  delay(100);
}

void receiveEvent(int howMany) {
  if (0 < Wire.available()) {   // receive a message
    int x = Wire.read();        // receive byte as an integer
    myservo.write(x);           // sets the servo position 
    delay(100);                 // waits for the servo to get there
    Serial.print(x); 
  }      
}
