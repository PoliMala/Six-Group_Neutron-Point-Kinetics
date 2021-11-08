#                     Neutron Point Kinetics simulation code   

Six Group of delayed neutron precursors implemented

It hence deal with a stiff problem due to the differences in evolution time scale between prompt and delayed neutron. For this purpose one of the best choice I discovered till now is the lsoda solver, which implement a dynamical choice of the solver method that adapts as best as possible to the model one wants to simulate.

Time index pointer version stands for the fact that, in view of coupling different models for the different physical processes one is interested in, the reactivity is managed by means an array which is loaded every time step.
Indeed a simpler version can be made just by means of functions properly defined but, in this way one can couple the different feedback processes just loading in the solving cycle the reactivity formula 
________________________________________________________________________________
#  MAIN Files Index

main.py   -> main program, performs the simulation and saves the csv with the results as well as the outcoming plots

KINmodule -> module containing the definition of the problem and the function to be implemented in order to solve the problem

KINnt     -> .png file with the outcoming neutron population time profile

KINct     -> .png file with the outcoming set of precursors group profile

.replit   -> config file for Repl usage - DO NOT DELETE

________________________________________________________________________________
#  SECONDARY Files Index

KINnt-promprjump -> .png file where displayed the rapid neutron jump due to prompt neutron
