%% Project 0: Simple Communication Link Simulation
% Alon S. Levin
% ECE-408: Wireless Communications
% Adapted from an assigment presented in ECE-300

clear, clc, close all;

%% Parameters
% Simulation parameters
numIter = 100;              % Number of iterations of the simulation
numSyms = 1000;             % Number of symbols per packet
SNR_Vec = 0:2:16;           % Vector containg various SNRs
lenSNR  = length(SNR_Vec);  % Number of SNRs we're testing at

% Modulation parameters
M = 16;                             % Modulation order
k = log2(M);                        % Number of bits per symbol
constellation = qammod(0:M-1,M);    % QAM reference constellation 

% Channel parameters
channel_ID = 1;                 % Choosing a channel mode, in increasing difficulty

% Equalizer parameters
order     = 12;                 % Order of the LMS equalizer
mu        = 0.7;                % Step-size of the LMS equalizer
equalizer = dsp.LMSFilter('Length', order, 'Method', 'Normalized LMS', 'StepSize', mu);
percTrain = 0.1;                % Percentage of bits transmitted that are to be trained on

% Trellis parameters
% 2/3 Feedforward Convolutional Encoder Trellis obtained from examples
% given by Mathworks
trellis          = poly2trellis([5 4], [23 35 0; 0 5 13]);
traceback_length = 20;          % Longer traceback yields better results, but longer computation time.

%% Performance Metrics
% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

%% Defining a Channel
% Note: Parameters work for cases 0-1, not yet for 2-3.
switch channel_ID
    case 0      % No channel
                % No ISI
        chan = 1;
    case 1      % Somewhat invertible channel impulse response
                % Moderate ISI
        chan = [1 .2 .4];
    case 2      % Worse channel
                % Severe ISI
        chan = [0.227 0.460 0.688 0.460 0.227];
    case 3      % Time-varying Rayleigh multipath channel
        ts = 1/1000;
        chan = rayleighchan(ts,1);
        chan.pathDelays    = [0 ts 2*ts];
        chan.AvgPathGaindB = [0 5 10];
        chan.StoreHistory  = 1;
    otherwise   % If improper ID, assume no channel
        chan = 1;
end

%% Simulation
% Run the simulation numIter amount of times

for ticker_iteration = 1:numIter
    %% Train Equalizer
    % Trains a new set of bits every iteration, converts to decimal and
    % then modulates.
    if channel_ID ~= 0  % Turn off equalizer when no channel
        % Generate training bits
    	train_bits = randi(2, [ceil(numSyms*percTrain), log2(M)]) - 1;
    	train_de   = bi2de(train_bits);
        
        % Modulate training bits
    	trainSig_de = qammod(train_de,M);
        
    	% Transmit training bits through channel
        if isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            rxSig_de = filter(chan,trainSig_de);
        else
            rxSig_de = filter(chan,1,trainSig_de);  % Apply the channel.
        end
        
        % Train equalizer
    	[~, ~, wts] = equalizer(rxSig_de, trainSig_de);
   	 
    end

    %% Generate a Message
    % New messages must be generated each iteration
    bits = randi(2, [ceil(numSyms*(1-percTrain)), log2(M)]) - 1;
    
    %% Apply Convolutional Coding, Prepare for Transmission
    code = convenc(bits(:), trellis);
    
    % Data needs to be reshaped into a matrix in order to allow for decimal
    % conversion and modulation.
    bin_msg = reshape(code, [length(code)/log2(M), log2(M)]);
    msg     = bi2de(bin_msg);
    tx      = qammod(msg, M);
    
    %% Simulate at each SNR Value
    % We send the message through the channel at different SNR values.
    for ticker_SNR = 1:lenSNR
        %% Transmit through Channel
        if channel_ID == 0
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan)                 % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
    
        %% Add AWGN
        % As per Prof. Keene's recommendation on this project, an SNR
        % correction should be added for higher-order modulation schemes
        SNR_correction = 10*log10(k);
        txNoisy = awgn(txChan, SNR_correction + SNR_Vec(ticker_SNR), 'measured');
        
        %% Apply Equalizer
        % Only if channel is present
        if channel_ID == 0
            eqSig = txNoisy;
        else
            eqSig = filter(wts, 1, txNoisy);
        end
        
        %% Demodulate the Signal
        % Here we also convert to binary column format to allow for Viterbi
        rx         = qamdemod(eqSig, M);
    	bin_rx     = de2bi(rx);
    	bin_rx_col = bin_rx(:);
        
    	decData = vitdec(bin_rx_col, trellis, traceback_length, 'trunc', 'hard');
    	rxMSG   = reshape(decData,[ceil(numSyms*(1-percTrain)), log2(M)]);
        
        %% Compute Bit Error
        [~, berVec(ticker_iteration, ticker_SNR)] = biterr(bits, rxMSG); 
    end
    
end

%% Results

% Constellation Diagram
h = scatterplot(txNoisy,1,0,'c.');
hold on
scatterplot(eqSig,1,0,'b.',h)
scatterplot(constellation,1,0,'r*',h);
legend('Received Signal','Equalized Signal','Constellation')
title(['LMS Equalizer, SNR = ', num2str(SNR_Vec(7)), ', Iteration ', num2str(ticker_iteration)])
hold off

% Compute and plot the mean BER
ber = mean(berVec,1);
figure
semilogy(SNR_Vec, ber)

% Compute the theoretical BER for this scenario
% THIS IS ONLY VALID FOR BPSK!
% YOU NEED TO CHANGE THE CALL TO BERAWGN FOR DIFF MOD TYPES
% Also note - there is no theoretical BER when you have a multipath channel
berTheory = berawgn(SNR_Vec,'qam',M);
hold on
semilogy(SNR_Vec,berTheory,'r')
legend('Empirical', 'Theoretical')
title('Bit Error Rate, Empirical vs Theoretical')
xlabel('SNR')
ylabel('Bit Error Rate')

bitRate = (1-percTrain)* size(bits(:)) / size(code);

fprintf("BER at 12dB = %d, BitRate = %d \n", ber(7), bitRate);