%% MIMO-OFDM Project
% Alon S. Levin
% ECE-408: Wireless Communications
% Spring 2020
% Part 2: OFDM

%% Prepare Environment
clear, clc, close all       % Clear all current variables and outputs
format compact              % Prepare command line output

%% Simulation Parameters

% System parameters
M = 4;                  % Modulation order: QAM
numSyms = 48e2;         % Number of transmitted symbols

% Channel parameters
numChan = 3;                                % Number of channels
numSNR = 20;                                % Number of SNR values to check
EbNo_vect = linspace(-10, numSNR, 20);      % Eb/No vector
SNR_vect = EbNo_vect + 10*log10(64/80); 	% SNR vector

% Rayleigh Channel parameters
Ts = 1e-3;                                  % Sample rate for Rayleigh Channel
Fd = 0;                                     % Maximum Doppler Shift
tau = [0 Ts/5 Ts/3 Ts];                     % Path delays
pdb = [0 -2 -3.33 -10];                     % Average path gain


%% Initialize BER Vectors
% BER_ZF   --> Zero-forcing
% BER_MMSE --> Minimum Mean Squared Error
BER_ZF   = nan(numChan, numSNR);
BER_MMSE = nan(numChan, numSNR);

%% Zero-Forcing
for channel_ticker = 1:numChan
    % Build the Rayleigh channel
    rayleighchan = comm.RayleighChannel(...
        'SampleRate', Ts, ...
        'PathDelays', tau, ...
        'AveragePathGains', pdb, ...
        'MaximumDopplerShift', Fd, ...
        'RandomStream','mt19937ar with seed', ...
        'Seed', randi(1e7) ...
    );

    % Generate a transmit signal
    tx_syms = randi([0, M-1], 1, numSyms);
    tx_mod = qammod(tx_syms, M);
    
    % Convert signal to frames
    tx_ofdm_frames = makeOFDMframe(tx_mod);
    numFrames = size(tx_ofdm_frames, 2);
    
    % Apply IFFT
    tx_postifft = ifft(tx_ofdm_frames, 64);
    
    % Add cyclic prefix
    tx_withcp = [tx_postifft(49:64,:); ...
                 tx_postifft];
    
    % Transmit through the Rayleigh channel
    tx_rayleigh = zeros(size(tx_withcp));
    channel_state = zeros(size(tx_withcp));
    for frame_ticker = 1:numFrames
        tx_rayleigh(:, frame_ticker) = rayleighchan(tx_withcp(:, frame_ticker));
        channel_state(:, frame_ticker) = rayleighchan(ones(80,1));
    end
    
    % Add AWGN per SNR values, post-code, find BER
    for SNR_ticker = 1:numSNR
        currSNR = SNR_vect(SNR_ticker);

        % Generate AWGN
        awgn = 10^(-currSNR/20) * sqrt(1/2) * ...
            (randn(80, numFrames) + 1j*randn(80, numFrames));

        % Add AWGN to the channel
        tx_awgnchannel = tx_rayleigh + awgn;
        
        % Strip cyclic prefix, apply FFT
        rx_nocp = tx_awgnchannel(17:end, :);
        rx_fft = fft(rx_nocp, 64);
        
        % Apply zero-forcing equalizer
        rx_zf = rx_fft ./ channel_state(17:end, :);
        
        % Restore data
        rx_data = rx_zf([...
            6:10, ...   % -26 --> -22
            12:24, ...  % -20 --> -08
            26:31, ...  % -06 --> -01
            33:38, ...  %  01 -->  06
            40:52, ...  %  08 -->  20
            54:58], :); %  22 -->  26
        rx_data_reshaped = reshape(rx_data, 1, []);
        
        % Demodulate, compute BER
        rx_demod = qamdemod(rx_data_reshaped, M);
        BER_ZF(channel_ticker, SNR_ticker) = mean(rx_demod ~= tx_syms);
    end
end

%% MMSE
for channel_ticker = 1:numChan
    % Generate a transmit signal
    tx_syms = randi([0, M-1], 1, numSyms);
    tx_mod = qammod(tx_syms, M);
    
    % Convert signal to frames
    tx_ofdm_frames = makeOFDMframe(tx_mod);
    numFrames = size(tx_ofdm_frames, 2);
    
    % Apply IFFT
    tx_postifft = ifft(tx_ofdm_frames, 64);
    
    % Add cyclic prefix
    tx_withcp = [tx_postifft(49:64,:); ...
                 tx_postifft];
    
    % Build the Rayleigh channel
    rayleighchan = comm.RayleighChannel(...
        'SampleRate', Ts, ...
        'PathDelays', tau, ...
        'AveragePathGains', pdb, ...
        'MaximumDopplerShift', Fd, ...
        'RandomStream','mt19937ar with seed', ...
        'Seed', randi(1e7) ...
    );
    
    % Transmit through the Rayleigh channel
    tx_rayleigh = zeros(size(tx_withcp));
    channel_state = zeros(size(tx_withcp));
    for frame_ticker = 1:numFrames
        tx_rayleigh(:, frame_ticker) = rayleighchan(tx_withcp(:, frame_ticker));
        channel_state(:, frame_ticker) = rayleighchan(ones(80,1));
    end
    
    % Add AWGN per SNR values, post-code, find BER
    for SNR_ticker = 1:numSNR
        currSNR = SNR_vect(SNR_ticker);

        % Generate AWGN
        awgn = 10^(-currSNR/20) * sqrt(1/2) * ...
            (randn(80, numFrames) + 1j*randn(80, numFrames));

        % Add AWGN to the channel
        tx_awgnchannel = tx_rayleigh + awgn;
        
        % Strip cyclic prefix, apply FFT
        rx_nocp = tx_awgnchannel(17:end, :);
        rx_fft = fft(rx_nocp, 64);
        
        % Apply zero-forcing equalizer
        norm = conj(channel_state(17:end,:)) .* channel_state(17:end,:) + 10^(-currSNR/20);
        rx_mmse = rx_fft .* conj(channel_state(17:end,:)) ./ norm;
        
        % Restore data
        rx_data = rx_mmse([...
            6:10, ...   % -26 --> -22
            12:24, ...  % -20 --> -08
            26:31, ...  % -06 --> -01
            33:38, ...  %  01 -->  06
            40:52, ...  %  08 -->  20
            54:58], :); %  22 -->  26
        rx_data_reshaped = reshape(rx_data, 1, []);
        
        % Demodulate, compute BER
        rx_demod = qamdemod(rx_data_reshaped, M);
        BER_MMSE(channel_ticker, SNR_ticker) = mean(rx_demod ~= tx_syms);
    end
end

%% Plot BERs
figure
semilogy(EbNo_vect, mean(BER_ZF)', '-.*', 'LineWidth', 2)
hold on
semilogy(EbNo_vect, mean(BER_MMSE)', '-.^', 'LineWidth', 2)
title('BER for QAM OFDM System')
grid on
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
legend({'Zero-Forcing', 'MMSE'}, 'Location', 'southwest')