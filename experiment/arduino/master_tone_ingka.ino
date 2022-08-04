//Código para el arduino maestro
//Se comunica con Matlab y hace la parte sonora del experimento
//Envía parámetros mecánicos al arduino esclavo
//Mueve al Servo entre dos bips


#include <Tone.h>
#include <stdlib.h>
#include <Wire.h>

#define STIMPIN 13	//stimulus output to pin 13
#define RESPPIN 12	//response feedback output to pin 12
#define INPUTPIN 11	//response interruptor between pin 4 and 5V
#define SENSORPIN 5
#define INGKAPIN 52


#define STIM_DURATION 50	//stimulus duration (milliseconds)
#define ANTIBOUNCE (0.5*isi)//minimum interval between responses (milliseconds)


//---------------------------------------------------------------------
//definition of global variables

//trial parameters
unsigned int isi,n_stim;
int perturb_size;
unsigned int perturb_bip;
float baseline_tau;

//general variables
int i,resp;
unsigned int stim_number,resp_number,force_number,asyn_number,delta_number,isi_number;
unsigned int *asyn_stimnum;
unsigned long t,prev_stim_t,prev_resp_t,prev_force_t;
char *event_name;
unsigned int *event_number,*event_coord;
boolean *event_used_flag;
int asynchrony,*asyn;
int delta_isi;
unsigned long *event_time,last_resp_time,last_stim_time;
unsigned int event_counter;
int last_resp_indx,last_force_indx,last_stim_indx,prevlast_stim_indx,prevlast_resp_indx,prevlast_force_indx;
int TH=50;                 // Valor del umbral del sensor
int INGKA = LOW;

//Servo variables
int servo_ini_pos,mechanical_size,mechanical_bip;
boolean mech_flag;

char message[20];
Tone stim_tone,resp_tone;
boolean allow,perturb_flag,perturb_mark_flag,asynchrony_defined,missing_resp,extra_resp,bad_asyn;


//---------------------------------------------------------------------
int memorialibre=0;
extern int __bss_end;
extern void *__brkval;

int get_free_memory() {
	int free_memory;

	if((int)__brkval == 0)
		free_memory = ((int)&free_memory) - ((int)&__bss_end);
	else
		free_memory = ((int)&free_memory) - ((int)__brkval);
	return free_memory;
}

//---------------------------------------------------------------------
//print a line avoiding "carriage return" of Serial.println()
void serial_print_string(char *string) {
	Serial.print(string);
	Serial.print("\n");
	return;
}

//---------------------------------------------------------------------
//decimal (8 bits) to binary conversion
//result: least significant bit (even/odd) to the right
void fast_d2b(unsigned char dec, unsigned char *c) {
	char i;
	for(i=0; i<8; i++)
		c[7-i] = (dec >> i) & 0x1;
	return;
}


//---------------------------------------------------------------------
//parse data from serial input
//input data format: eg "I500;N30;P-10;B15;E5;X"
void parse_data(char *line) {
	char field[10];
	int n,data;
	//scan input until next ';' (field separator)
	while (sscanf(line,"%[^;]%n",field,&n) == 1) {
		data = atoi(field+1);
		//parse data according to field header
		switch (field[0]) {
			case 'I':
				baseline_tau = data;
				break;
			case 'N':
				n_stim = data;
				break;
			case 'p':
				perturb_bip = data;
				break;
			case 'P':
				perturb_size = data;
				break;
            case 'S':
				servo_ini_pos = data;
				break;
            case 'M':
				mechanical_size = data;
				break;
            case 'm':
				mechanical_bip = data;
				break;
			default:
				break;
		}
		line += n;
		if (*line != ';')
			break;
		while (*line == ';')
			line++;
	}
	return;
}
//---------------------------------------------------------------------
void get_parameters() {
	char line[45],i,aux='0';
	i = 0;

	//directly read next available character from buffer
	//if flow ever gets here, then next available character should be 'I'
	aux = Serial.read();

	//read buffer until getting an X (end of message)
	while (aux != 'X') {
		//keep reading if input buffer is empty
		while (Serial.available() < 1) {}
		line[i] = aux;
		i++;
		aux = Serial.read();
	}
	line[i] = '\0';					//terminate the string

	//just in case, clear incoming buffer once read
	Serial.flush();
	//parse input chain into parameters
	parse_data(line);
	return;
}

//---------------------------------------------------------------------
//get keyword before reading parameters
void get_keyword() {
	boolean allow = false;
	char keywrd[5] = "EEEE";

	//read input buffer until keyword "ARDU;"
	while (allow == false) {
		while (Serial.available() < 1) {
			//allow user to tap while waiting for data from computer
			resp = analogRead(SENSORPIN);
			if (resp >= TH)
				resp_tone.play(NOTE_D5,STIM_DURATION);
		}
		//read input buffer one at a time
		if (keywrd[0] == 'A' && keywrd[1] == 'R'
			&& keywrd[2] == 'D' && keywrd[3] == 'U' && keywrd[4] == ';') {
			//only combination allowed, in this specific order
			allow = true;
		}
		else {
			//move buffer one step up
			keywrd[0] = keywrd[1];
			keywrd[1] = keywrd[2];
			keywrd[2] = keywrd[3];
			keywrd[3] = keywrd[4];
			keywrd[4] = Serial.read();
		}
	}
	return;
}

//---------------------------------------------------------------------
void setup() {
  Serial.begin(9600);     //USB communication with computer
  Wire.begin(); // join i2c bus (address optional for master)
  stim_tone.begin(STIMPIN);   //stimulus output
  resp_tone.begin(RESPPIN);   //feedback output
  pinMode(INPUTPIN,INPUT);    //input pin pulled-down to GND by actual resistor
  pinMode(INGKAPIN, OUTPUT);      // sets the digital pin as output
  digitalWrite(INGKAPIN,LOW);    //just in case
  digitalWrite(INPUTPIN,LOW);   //just in case
  allow = false;


}
//---------------------------------------------------------------------
//main loop
void loop() {
	if (allow == false) {
		//constantly check buffer for keyword
		get_keyword();
		//ok, keyword found
		allow = true;
		missing_resp = false;
		extra_resp = false;
		bad_asyn = false;
		get_parameters();

		isi = baseline_tau;
		delta_isi = 0;
		prev_stim_t = 0;
		prev_resp_t = 0;
		prev_force_t = 0;
		stim_number = 0;
		resp_number = 0;
		asyn_number = 0;
		force_number = 0;
		event_counter = 0;
		last_resp_indx = -1;
		last_force_indx = -1;
		last_stim_indx = -1;
		prevlast_stim_indx = -1;
		prevlast_resp_indx = -1;
		prevlast_force_indx = -1;
		perturb_flag = false;
		perturb_mark_flag = false;
		mech_flag=false;

		event_name = (char*) calloc(10*n_stim,sizeof(char)); //labels para estimulos, respuestas, asincronias y delta_isis
		event_number = (unsigned int*) calloc(10*n_stim,sizeof(unsigned int)); //numera cada tipo de evento independientemente
		event_coord = (unsigned int*) calloc(10*n_stim,sizeof(unsigned int)); //alinea respuestas y asincronias al estimulo correspondiente
		event_time = (unsigned long*) calloc(10*n_stim,sizeof(unsigned long)); //valor de cada evento (tiempos para stim y resp, duraciones para asyn y delta_isi)
		event_used_flag = (boolean*) calloc(10*n_stim,sizeof(boolean)); //marca los eventos stim/resp que ya fueron usados para calcular una asincronia

		for (i=0; i<10*n_stim; i++) {
			event_used_flag[i] = false;
		}
   
   //reset Servo position
   // Wire.beginTransmission(8);    // transmit to device #8
   // Wire.write(servo_ini_pos);    // send the signal to move Servo
   // Wire.endTransmission();       // stop transmitting
   // delay(500);		
	}
	//start trial
	else {		
		t = millis();

		//turn on stimulus
		if ((t - prev_stim_t) >= isi && stim_number < n_stim) {
			//if last_stim_indx not used and previous used then break (i.e. missing response)
			if ((stim_number>=2 && event_used_flag[last_stim_indx] == false && event_used_flag[prevlast_stim_indx] == true) || (stim_number==5 && prev_resp_t==0 )){	
				missing_resp = true;
				stim_tone.play(NOTE_C6,STIM_DURATION);
			}
			else {
				stim_number++;
				stim_tone.play(NOTE_A4,STIM_DURATION);
				// señal al ingka
				INGKA = HIGH;
				digitalWrite(INGKAPIN,INGKA);
        		
				//store event data
				event_name[event_counter] = 'S';
				event_number[event_counter] = stim_number;
				event_coord[event_counter] = stim_number;
				event_time[event_counter] = t;
				prevlast_stim_indx = last_stim_indx;
				last_stim_indx = event_counter;
				event_counter++;
				prev_stim_t = t;			
						}
		}
		
		if ((t - prev_stim_t) >= STIM_DURATION  && INGKA == HIGH) {
              //INGKA sync pulse OFF
              INGKA =LOW;
              digitalWrite(INGKAPIN,INGKA);
        }	
		
		//spatial perturbation
		if (stim_number==mechanical_bip && (t - prev_stim_t) >= 0.5*isi && mech_flag==false ){
			mech_flag=true;
			Wire.beginTransmission(8); // transmit to device #8
			Wire.write(mechanical_size); // send the signal to move Servo
			Wire.endTransmission();    // stop transmitting
		}
		
		//read response
		if ((t - prev_resp_t) > ANTIBOUNCE) {
			resp = analogRead(SENSORPIN);
			if (resp >= TH){
				//if last_resp_indx not used and previous used then break (i.e. extra response)
				if (resp_number>=2 && event_used_flag[last_resp_indx] == false && event_used_flag[prevlast_resp_indx] == true) {
					extra_resp = true;
					stim_tone.play(NOTE_C6,STIM_DURATION);
				}
				else {
					resp_number++;
					force_number++;
					resp_tone.play(NOTE_D5,STIM_DURATION);
					
					//store event data
					event_name[event_counter] = 'R';
					event_number[event_counter] = resp_number;
					event_time[event_counter] = t;
					prevlast_resp_indx = last_resp_indx;
					last_resp_indx = event_counter;
					event_counter++;
					prev_resp_t = t;
					
					//store event data
					event_name[event_counter] = 'F';
					event_number[event_counter] = force_number;
					event_time[event_counter] = resp;
					//prevlast_force_indx = last_force_indx;
					//last_force_indx = event_counter;
					        event_coord[event_counter] = event_number[last_stim_indx];
        event_coord[last_resp_indx] = event_number[last_stim_indx];
        event_counter++;
					prev_force_t = t;
				}
			}
		}
		
		//only test if the previous two events were not used already
		//asynchrony defined only if a pair R-S or pair S-R
		if ((last_stim_indx >= 0 && last_resp_indx >= 0) && event_used_flag[last_resp_indx] == false && event_used_flag[last_stim_indx] == false) {
			last_resp_time = event_time[last_resp_indx];
			last_stim_time = event_time[last_stim_indx];
			asynchrony = last_resp_time - last_stim_time;

			if (resp_number == 1 && abs(asynchrony) >= 0.5*isi) {
				//just skip and wait for next stimulus
			}
			else if (resp_number > 1 && abs(asynchrony) >= 0.5*isi) {
				bad_asyn = true;
				stim_tone.play(NOTE_C6,STIM_DURATION);
				asyn_number++;
				event_used_flag[last_resp_indx] = true;
				event_used_flag[last_stim_indx] = true;
				event_name[event_counter] = 'A';
				event_time[event_counter] = asynchrony;
				event_number[event_counter] = asyn_number;
				event_coord[event_counter] = event_number[last_stim_indx];
				event_coord[last_resp_indx] = event_number[last_stim_indx];
				event_counter++;
			}
			else {
				//first R-S or S-R pair ok
				asyn_number++;
				event_used_flag[last_resp_indx] = true;
				event_used_flag[last_stim_indx] = true;
				event_name[event_counter] = 'A';
				event_time[event_counter] = asynchrony;
				event_number[event_counter] = asyn_number;
				event_coord[event_counter] = event_number[last_stim_indx];
				event_coord[last_resp_indx] = event_number[last_stim_indx];
				event_counter++;

				//step change
				if (stim_number >= perturb_bip) {
					if (stim_number == perturb_bip) {
						delta_isi = (int) perturb_size;
					}
					else {
						delta_isi = 0;
					}
				} 			
				//change isi
				isi = isi + delta_isi;			
			}
		}
		//end trial
		//allow one more period (without stimulus)
		if ((stim_number == n_stim && (t - prev_stim_t) >= isi) || missing_resp == true || extra_resp == true || bad_asyn == true) {

      Serial.println("A");
      
			for (i=0; i<event_counter; i++) {
      
				//print events
				sprintf(message,"%c %d(%d): %ld;",
					event_name[i],event_number[i],event_coord[i],event_time[i]);
				serial_print_string(message);
			}
			
			//print error message
			if (missing_resp == true){
				sprintf(message,"M");
				serial_print_string(message);
			}
			if (extra_resp == true){
				sprintf(message,"X");
				serial_print_string(message);
			}
			if (bad_asyn == true){
				sprintf(message,"B");
				serial_print_string(message);
			}
			
			//send END of message
			Serial.println("E");	
			allow = false;
			free(event_name);
			free(event_number);
			free(event_time);
			free(event_coord);
			free(event_used_flag);
			
			delay(500);
			Wire.beginTransmission(8);    // transmit to device #8
			Wire.write(servo_ini_pos);    // send the signal to move Servo
			Wire.endTransmission();       // stop transmitting			
		}
	}
}
