# EletrophysiologyAnalysis
Code for population analysis of electrophysiological data

Eletrophysiological code used for population analysis of acute superior collicular recordings under subcutaneous medetomidine sedation. 


Electrophysiology probe used: 
64 channel (8 shanks × 8 sites) silicon probe (Buzsaki64, NeuroNexus, Ann Arbor, MI). 

Channel map order:
orderVector= [48,35,44,46,43,38,40,39,42,54,33,34,37,52,50,36,49,63,56,55,51,61,58,53,57,41,59,47,64,45,62,60,23,7,17,5,19,2,6,4,1,15,9,10,3,13,11,8,12,24,32,31,14,27,30,16,29,18,20,22,28,21,25,26]

Acquisition software: Open Ephys

Each run consists on 10 stimulation cycles.
One cycle: Stimulation period (15 seconds) + Rest period (45 seconds)

The code calculates the power spectral density and power spectrum of each recording.
A section for correlations with other data type (in this case was fMRI data) is available.
The code for the calculation of a neurometirc curve based on power spectrum peaks is available.
