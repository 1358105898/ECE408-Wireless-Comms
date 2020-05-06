%% MIMO-OFDM Project
% Alon S. Levin
% ECE-408: Wireless Communications
% Spring 2020
% Part 3: MIMO-OFDM

%% Prepare Environment
clear, clc, close all       % Clear all current variables and outputs
format compact              % Prepare command line output

%% Simulation Parameters

% System parameters
M = 4;                  % Modulation order: QAM
numSyms = 48e2;         % Number of transmitted symbols
numTx = 2;              % Number of transmitters
numRx = 2;              % Number of receivers

% Channel parameters
numChan = 3;                                % Number of channels
numSNR = 20;                                % Number of SNR values to check
EbNo_vect = linspace(-10, numSNR, 20);      % Eb/No vector
SNR_vect = EbNo_vect + 10*log10(log2(M)); 	% SNR vector

% Rayleigh Channel parameters
Ts = 1e-3;                                  % Sample rate for Rayleigh Channel
Fd = 0;                                     % Maximum Doppler Shift
tau = [0 Ts/5 Ts/3 Ts];                     % Path delays
pdb = [0 -2 -3.33 -10];                     % Average path gain