clear; close all; clc

%% Simscape Controlling Script
% Run the following commands to simulate the 'helmoltz_tree_SSC.slx' file 
% and generate the ground truth signal as a MATLAB workspace variable.

% Single Resonator Parameters
V  = 0.10;    % m^3
l  = 10.0;    % m
S  = 100.0;   % m^2

rho = 1.2;    % kg/m^3
c   = 343.2;  % m/s

% Electrical Equivalent Model Parameters
R = (rho * c) / S; 
L = (rho * l) / S;
C = V / (rho * c^2);

% Sampling Frequency
fs = 96e3;

% Sampling Period
Ts = 1/fs;

% Simulation Duration
stop_time = 5;  % [seconds]

% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1).';

vin = zeros(length(t), 1);
vin(2) = 1;

% Define Simscape Input Signal (ssc_input)
ssc_input = [t, vin];

%% Model Simulation

model_name = "helmholtz_tree_SSC";

sim(model_name);