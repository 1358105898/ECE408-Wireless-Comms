function [BER_vector] = estimateBER_Alamouti(numRx, numBits, SNR_vector, numIter, f_m)
%estimateBER_Alamouti Runs an Alamouti simulation through a Rayleigh
%channel, and reports the BER_vector at each given SNR

% Prepare BER vector to store results
numSNR = length(SNR_vector);
BER_vector = zeros(numSNR, 1);

% Begin iterating
for iter_ticker = 1:numIter
    for SNR_ticker = 1:numSNR
        currSNR = SNR_vector(SNR_ticker);
        % Generate data, convert to BPSK, stored as column vector
        data_raw = randi([0, 1], numBits, 1);
        data_tx_linear = 2*data_raw - 1;
        
        % Split Tx data across two antennas, as per Alamouti's algorithm
        data_tx_split = reshape(data_tx_linear, [2, numBits/2]);
        data_tx_0 = [data_tx_split(1,:); ...
                    -data_tx_split(2,:)];
        data_tx_1 = [data_tx_split(2,:); ...
                     data_tx_split(1,:)];
        data_tx   = [data_tx_0(:), data_tx_1(:)];
        
        % Generate the channels, stored as column vectors
        chans = nan(numBits, 2*numRx);
        for channel_ticker = 1:(2*numRx)
            temp_chan = generateRayleighChannel(numBits/2, f_m).';
            chan_split = [temp_chan; temp_chan];
            chans(:, channel_ticker) = chan_split(:);
        end
        
        % Pass the data through the channels
        % Note the 3dB penalty for transmitting at half power
        data_rx = nan(numBits, numRx);
        for Rx_ticker = 1:numRx
            data_rx(:, Rx_ticker) = ...
                awgn(sum(data_tx .* chans(:, Rx_ticker:Rx_ticker+1), 2), ...
                     currSNR-3, 'measured');
        end
        
        % Split incoming stream into r0 and r1 vectors
        r0 = data_rx(1:2:numBits, :);
        r1 = data_rx(2:2:numBits, :);
        
        % Estimate symbol by multiplying channel conjugate and received
        % signal on each antenna
        s_tilde_0 = zeros(numBits/2, 1);
        s_tilde_1 = zeros(numBits/2, 1);
        for Rx_ticker = 1:numRx
            s_tilde_0 = s_tilde_0 + ...
                        r0(:, Rx_ticker) .* conj(chans(1:2:numBits, Rx_ticker)) + ...
                        conj(r1(:, Rx_ticker)) .* chans(1:2:numBits, Rx_ticker+1);
            s_tilde_1 = s_tilde_1 + ...
                        r0(:, Rx_ticker) .* conj(chans(1:2:numBits, Rx_ticker+1)) - ...
                        conj(r1(:, Rx_ticker)) .* chans(1:2:numBits, Rx_ticker);
        end
        
        % Merge the symbols together
        s_tilde_split = [s_tilde_0.'; ...
                         s_tilde_1.'];
        s_tilde = s_tilde_split(:);
        
        % Estimate the symbol, computer BER
        data_est = (s_tilde > 0);
        BER_vector(SNR_ticker) = BER_vector(SNR_ticker) + mean(xor(data_raw, data_est));
    end
end

% Normalize BER vector
BER_vector = BER_vector / numIter;

end


