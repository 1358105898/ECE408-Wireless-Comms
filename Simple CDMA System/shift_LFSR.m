function [newLFSR, output] = shift_LFSR(LFSR, taps, numShifts)
%shift_LFSR shifts the LFSR using the taps listed, and returns the new LFSR
%state and the output (equivalent to LFSR(1))
% Note: assume taps are defined properly for the LFSR type

% Pre-allocate outputs
output = nan(1, numShifts);
newLFSR = LFSR;

% Length of LFSR
numStages = length(LFSR);

% Define the positions of the taps
taps_positions = ismember(1:numStages, taps);

% Begin shifting
for ticker = 1:numShifts
    % Obtain the new entry
    new_entry = mod(sum(taps_positions.*LFSR), 2);

    % Shift
    newLFSR(1:numStages-1) = LFSR(2:end);
    newLFSR(numStages) = new_entry;
    output(ticker) = newLFSR(1);
    
    LFSR = newLFSR;
end

end

