% default parameters
runID = 'byoplanet';   % change run identifier for every parameter test
nop   =  100;          % plot/print output figures every nop time steps

seed  =  15;           % set the random seed for reproducible results

tend  =  1000;         % stopping time for testing, higher for full run

N     =  1000;         % initial number of planetesimals

MPls  =  0.03;         % mean planetesimal mass [Earth Masses]
MGgt  =  300;          % mass of gas giant [Earth Masses]
MStr  =  3e5;          % mass of central star [Earth Masses]

cls   =  0.03;         % effective collision radius [AU]