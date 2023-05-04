# channel flows basics
This repository contains codes for calculating some basic properites for channel flow using one instantaneous velocity profile. For time-averaged properities, you should attain sufficient amount of snapshots. 
File 'velocity_data.mat' contains one instantaneous velocity profile, u is the streamwise velocity, v is the spanwise velocity, w is the wall-normal velocity. The mean flow and its first order derivative with respect to the wall-normal direction are also given.

Code 'calculate_PDN_F.m' calculates the production, dissipation and nonlinear transfer for the entire channel flow. The values at specific wall-normal heights can also be calculated by slightly modifying the codes.

Code 'calculate_M.m' calculates the source-recipient nonlinear transfer between two scales. 

Code 'calculate_swirling_strength.m' calculates the swirling strength and implements the visualisation. 
