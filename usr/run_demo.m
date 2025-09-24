% prepare workspace
clear; close all; clc;

% load default parameters
par_default;

runID = 'byoplanet';   % change run identifier for every parameter test
nop   =  100;          % plot/print output figures every nop time steps

seed  =  15;           % set the random seed for reproducible results

tend  =  1000;         % stopping time for testing, higher for full run

N     =  10000;        % initial number of planetesimals

MPls  =  0.01;         % mean planetesimal mass [Earth Masses]
MGgt  =  300;          % mass of gas giant [Earth Masses]
MStr  =  3e5;          % mass of central star [Earth Masses]

cls   =  0.01;         % effective collision radius [AU]

% run model
run('../src/main');