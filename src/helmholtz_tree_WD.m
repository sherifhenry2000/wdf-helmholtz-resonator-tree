%% Sherif Henry Mahroos Henry
clear; close all; clc

%% Simscape File for Ground Truth
load('ssc_output.mat')

%% Sampling Frequency
fs = 96e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 5;  % [seconds]

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1).';

%% Input Signal (Kronecker Delta)

vin = zeros(length(t), 1);
vin(2) = 1;

%% Input Signal (Noise Excitation)
% 
% % Reproducible random seed
% rng(42);   
% 
% % Windowed noise (20 ms)
% dur = 20e-3;
% Nw = round(dur * fs);
% 
% vin = zeros(length(t), 1);
% 
% vin(1:Nw) = randn(Nw, 1) .* hann(Nw);
% vin = vin ./ max(abs(vin));

%% Circuit Parameters
% Helmholtz Resonator Tree

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

[R00, R10, R11, R20, R21] = deal(R);
[L00, L10, L11, L20, L21] = deal(L);
[C00, C10, C11, C20, C21] = deal(C);

%% Setting of Free Parameters (Adaptation Conditions)
Z2 = R00;
Z16 = R10; 
Z10 = R11; 
Z28 = R20;
Z22 = R21;

Z4 = 2*L00/Ts;
Z14 = 2*L10/Ts; 
Z12 = 2*L11/Ts; 
Z26 = 2*L20/Ts; 
Z24 = 2*L21/Ts; 

Z6 = Ts/(2*C00);
Z19 = Ts/(2*C10);
Z11 = Ts/(2*C11);
Z27 = Ts/(2*C20);
Z23 = Ts/(2*C21);

Z25 = Z28 + Z26 + Z27;
Z21 = Z22 + Z24 + Z23;
Z9  = Z10 + Z12 + Z11;

Z18 = Z25; Z20 = Z21; 
Z17 = 1 / ( (1/Z18) + (1/Z19) + (1/Z20) );

Z15 = Z17; 
Z13 = Z16 + Z14 + Z15;

Z8 = Z9; Z7 = Z13; 
Z5 = 1 / ( (1/Z6) + (1/Z7) + (1/Z8) );

Z3 = Z5; 
Z1 = Z2 + Z4 + Z3;


%% Computing Scattering Matrices
% Followed these rules:
% For Series:   S = I - alpha * [1]^T (where alpha_j = 2*Z_j / sum(Z))
% For Parallel: S = [1] * alpha^T - I (where alpha_j = 2*G_j / sum(G))

% Adaptor 25 (4-port Series)
Z_A25 = [Z25, Z28, Z26, Z27];
alpha_A25 = 2 * Z_A25 / sum(Z_A25);
S_A25 = eye(4) - alpha_A25.' * ones(1,4); % FIXED TRANSPOSE

% Adaptor 21 (4-port Series)
Z_A21 = [Z21, Z22, Z24, Z23];
alpha_A21 = 2 * Z_A21 / sum(Z_A21);
S_A21 = eye(4) - alpha_A21.' * ones(1,4); 

% Adaptor 9 (4-port Series)
Z_A9 = [Z9, Z10, Z12, Z11];
alpha_A9 = 2 * Z_A9 / sum(Z_A9);
S_A9 = eye(4) - alpha_A9.' * ones(1,4); 

% Adaptor 17 (4-port Parallel) 
G_A17 = [1/Z17, 1/Z18, 1/Z19, 1/Z20];
alpha_A17 = 2 * G_A17 / sum(G_A17);
S_A17 = ones(4,1) * alpha_A17 - eye(4); 

% Adaptor 13 (4-port Series)
Z_A13 = [Z13, Z16, Z14, Z15];
alpha_A13 = 2 * Z_A13 / sum(Z_A13);
S_A13 = eye(4) - alpha_A13.' * ones(1,4); 

% Adaptor 5 (4-port Parallel) 
G_A5 = [1/Z5, 1/Z6, 1/Z7, 1/Z8];
alpha_A5 = 2 * G_A5 / sum(G_A5);
S_A5 = ones(4,1) * alpha_A5 - eye(4); 

% Adaptor 1 (4-port Series)
Z_A1 = [Z1, Z2, Z4, Z3];
alpha_A1 = 2 * Z_A1 / sum(Z_A1);
S_A1 = eye(4) - alpha_A1.' * ones(1,4); 


%% Initialization of Waves
% Delay states for Capacitors and Inductors
a4_del = 0; a6_del = 0;   
a14_del = 0; a19_del = 0; 
a12_del = 0; a11_del = 0; 
a26_del = 0; a27_del = 0; 
a24_del = 0; a23_del = 0;


%% Initialization of Output Signals

iout = zeros(1, length(t));

%% Simulation Algorithm

tic
for n = 1 : length(t)

    % Forward Scan
    % Read Leaves
    a28 = 0; a26 = -a26_del; a27 = a27_del;
    a22 = 0; a24 = -a24_del; a23 = a23_del;
    a10 = 0; a12 = -a12_del; a11 = a11_del;
    a16 = 0; a14 = -a14_del; a19 = a19_del;
    a2  = 0; a4  = -a4_del;  a6  = a6_del;

    % Matrix Multiply: v_b = S * v_a 
    
    v_b_A25 = S_A25 * [0; a28; a26; a27];
    b25 = v_b_A25(1);
    
    v_b_A21 = S_A21 * [0; a22; a24; a23];
    b21 = v_b_A21(1);
    
    v_b_A9 = S_A9 * [0; a10; a12; a11];
    b9 = v_b_A9(1);

    a18 = b25; a20 = b21; 
    v_b_A17 = S_A17 * [0; a18; a19; a20];
    b17 = v_b_A17(1);

    a15 = b17;
    v_b_A13 = S_A13 * [0; a16; a14; a15];
    b13 = v_b_A13(1);

    a8 = b9; a7 = b13;
    v_b_A5 = S_A5 * [0; a6; a7; a8];
    b5 = v_b_A5(1);

    a3 = b5;
    v_b_A1 = S_A1 * [0; a2; a4; a3];
    b1 = v_b_A1(1);

    
    % Local Root Scattering

    a_source = b1;
    b_source = 2*vin(n) - a_source;
    % Backward Scan

    a1 = b_source;

    v_b_A1 = S_A1 * [a1; a2; a4; a3];
    b2 = v_b_A1(2); b4 = v_b_A1(3); b3 = v_b_A1(4);

    a5 = b3;
    v_b_A5 = S_A5 * [a5; a6; a7; a8];
    b6 = v_b_A5(2); b7 = v_b_A5(3); b8 = v_b_A5(4);

    a13 = b7;
    v_b_A13 = S_A13 * [a13; a16; a14; a15];
    b16 = v_b_A13(2); b14 = v_b_A13(3); b15 = v_b_A13(4);
    
    a9 = b8;
    v_b_A9 = S_A9 * [a9; a10; a12; a11];
    b10 = v_b_A9(2); b12 = v_b_A9(3); b11 = v_b_A9(4);

    a17 = b15;
    v_b_A17 = S_A17 * [a17; a18; a19; a20];
    b18 = v_b_A17(2); b19 = v_b_A17(3); b20 = v_b_A17(4);

    a25 = b18;
    v_b_A25 = S_A25 * [a25; a28; a26; a27];
    b28 = v_b_A25(2); b26 = v_b_A25(3); b27 = v_b_A25(4);

    a21 = b20;
    v_b_A21 = S_A21 * [a21; a22; a24; a23];
    b22 = v_b_A21(2); b24 = v_b_A21(3); b23 = v_b_A21(4);

    % Delay memory updates for the next sample
    a4_del = b4; a6_del = b6;
    a14_del = b14; a19_del = b19;
    a12_del = b12; a11_del = b11;
    a26_del = b26; a27_del = b27;
    a24_del = b24; a23_del = b23;
    % Read Output
    iout(n) = (a2 - b2) / (2 * Z2);

end
toc

%% Output Plots

% Predicted signal
iout = iout.';

% Obtain Ground Truth signal
iout_gt = gt(2, 1:end-1).';

% Time Domain Plots
figure
set(gcf, 'Color', 'w');
plot(t, iout, 'b', 'Linewidth', 2);
hold on
plot(t, iout_gt, 'r--', 'Linewidth', 2);
grid on;
xlabel('Time [s]','Fontsize',16,'interpreter','latex');
ylabel('Volume Velocity $[\mathrm{m}^3/\mathrm{s}]$','Fontsize',16,...
    'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Impulse Response','Fontsize',18,'interpreter','latex');

% Frequency domain Plots
nfft = length(t);
res = fs/nfft;
f = (0:nfft/2-1) * res;

tf    = fft(iout, nfft);
tf_gt = fft(iout_gt, nfft);

abs_tf    = abs(tf(1:nfft/2));
abs_tf_gt = abs(tf_gt(1:nfft/2));

figure
set(gcf, 'Color', 'w');
plot(f, 20*log10(abs_tf), 'b', 'Linewidth', 2);
hold on
plot(f, 20*log10(abs_tf_gt), 'r--', 'Linewidth', 2);
grid on;
xlim([20, 2000])
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('Magnitude $[\mathrm{dB}]$','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Magnitude Response','Fontsize',16,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, iout - iout_gt, 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}\,[\mathrm{m}^3/\mathrm{s}]$',...
    'Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Normalized Mean Squared Error (NMSE)

nmse = sum((iout - iout_gt).^2) / sum(iout_gt.^2);

fprintf('NMSE = %.6g\n', nmse);

%% Write Audio File
% 
% y = iout - mean(iout);
% y = 0.95 * y ./ max(abs(y));
% audiowrite('helmholtz_tree_noise_response.wav', y, fs);
% 
% % Listen to the audio file
% soundsc(y, fs)