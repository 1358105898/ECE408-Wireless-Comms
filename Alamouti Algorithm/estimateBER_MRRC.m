function [BER_vector] = estimateBER_MRRC(numRx, numBits, SNR_vector, numIter, f_m)
%estimateBER_MRRC Runs an MRRC simulation through a Rayleigh channel, and
%reports the BER_vector at each given SNR

% Prepare BER vector to store results
numSNR = length(SNR_vector);
BER_vector = zeros(numSNR, 1);

% Begin iterating
for iter_ticker = 1:numIter
    for SNR_ticker = 1:numSNR
        currSNR = SNR_vector(SNR_ticker);
        % Generate data, convert to BPSK, stored as column vector
        data_raw = randi([0, 1], numBits, 1);
        data_tx = 2*data_raw - 1;

        % Generate the channels, stored as column vectors
        chans = nan(numBits, numRx);
        for channel_ticker = 1:numRx
            chans(:, channel_ticker) = generateRayleighChannel(numBits, f_m);
        end
        
        % Pass the data through the channels
        data_rx = awgn(data_tx .* chans, currSNR, 'measured');
        
        % Estimate symbol by multiplying channel conjugate and received
        % signal on each antenna
        s_tilde = sum(data_rx .* conj(chans), 2);
        
        % Estimate the symbol, computer BER
        data_est = (s_tilde > 0);
        BER_vector(SNR_ticker) = BER_vector(SNR_ticker) + mean(xor(data_raw, data_est));
    end
end

% Normalize BER vector
BER_vector = BER_vector / numIter;

end

