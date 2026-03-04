# Five-Tank-Cascade-System

The porpouse of this project was to control five tank in cascade with multiple controllers in different state-feedback control structures:  centralized, decentralized and two different distributed schemes: the first is a string with uniderectional communication and the second is a string with bidirectional communication.
The system has to be considered discrete.

The controllers used are: simple LMI, LMI with placing eigenvalues in a Disk(alpha, rho_target), H2 control and Hinf control.

There are two folders: MATLAB and images

MATLAB:

IrrigationChannel: is the file where all the initializzation take place. It's necessary to run only this file to decide see the effect of the various controller on the system. By typing a value is possible to decide which controller we want to analyse.

LMI_DT_DeDicont: is the file where there is the LMI for this discrete systems

LMI_DT_Didk_Struct: is the file where is build the LMI with placing eigenvalues in a Disk(alpha, rho_target)

LMI_DT_H2: is the file where H2 controller is set up

LMI_DT_Hinf: is the file wheere Hinf controller is set up

circle: is to build a cirlce necessary to analyse where the eigenvalues are respect the unit circle

di_fixed_modes: computes the fixed modes of the system

simulate_closed_loop: the name says itself

images:

control_action: contains all the images of the system responses, each file is called nameofcontroller_h#, # is the number of level

eigenvalues: contains images of the the eigenvalues of each structure controlled inside the unit circle after each controller. So the files are called nameofcontroller_structurecontrolled

