Three main codes wer wirrten in order to analyze the Ising Model.

ising.cpp: contains the principal routinesin order to modelate the evolution of the system under periodic or free consitions for a "square grid". The program is easily run by compiling first with make and the running ising_temp.out L number_of_steps 1==Perdiodic_contitions/0==Free_considions compiling first with make and the running ising_temp.out number_of_steps L 1==Save_energy/0==Save_Magnetization

On the other hand the program ising_system.cpp was done to analzize some properties of markov chains in the ising model. Do make ising5 and the run as ising_system.out L number_of_steps

In order for the program t run create the directories Periodic_data and Free_data.

Now in order to plot and analyze the data the program Analysis.py was done in python, the different analysis made are commentented and are easily understanandable.
