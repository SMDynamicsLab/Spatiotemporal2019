# Matlab files

To perform the experiments mentioned in the paper run **experiment_manager_simetrico.m** with the desired parameters

This script requires:

- **SYS_id.m**: detects the Operating System.
- **loop_central_simetrico.m**: establishes communication with the *master* Arduino, send it parameters, and receives data.
- **parameters_calibration_exp.m**: takes the data from the calibration stage and performs a linear fit to get the height required for a certain temporal perturbation.
- **hei2ang_servo.m**: translates the need height in centimeters to an angle in degrees for the platform position for the spatial perturbations.
- **save_bad_data.m**: save the data of bad trials.
- **save_good_data.m**: save the data of good trials.
- **f_assign_condition.m**: assign the condition number of the trial based on the types of perturbations performed.

# Arduino files

The folder *arduino* contains two files:

- **master_tone_ingka.ino** (for the *master* Arduino): receive the trial parameters (number of bips, size of perturbations. etc), controls the slave Arduino, perform the trial, and collect the data.
- **slave_servo.ino** (for the *slave* Arduino): moves the servo motor to the position indicated by the master Arduino.

More details about the setup available at:
*López, Sabrina Laura. (2021). Sincronización sensomotora. Perturbaciones temporales y espaciales en una tarea de finger tapping . (Tesis Doctoral. Universidad de Buenos Aires. Facultad de Ciencias Exactas y Naturales.)*
[link](https://bibliotecadigital.exactas.uba.ar/collection/tesis/document/tesis_n7018_Lopez)
