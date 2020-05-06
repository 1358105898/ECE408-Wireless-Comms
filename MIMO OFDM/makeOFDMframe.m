function [frame] = makeOFDMframe(data)
%makeOFDMframe uses the 802.11a frame parameters to turn generic qammod
%data into ofdm frame. In part based on the OFDM midterm project.

data_reshaped = reshape(data, 48, []);
numFrames = size(data_reshaped, 2);

frame = [zeros(5, numFrames); ...
         data_reshaped(1:5, :); ...    % -26 --> -22
         ones(1, numFrames); ...       % -21
         data_reshaped(6:18, :); ...   % -20 --> -08
         ones(1, numFrames); ...       % -07
         data_reshaped(19:24, :); ...  % -06 --> -01
         zeros(1, numFrames); ...
         data_reshaped(25:30, :); ...  %  01 -->  06
         ones(1, numFrames); ...       %  07
         data_reshaped(31:43, :); ...  %  08 -->  20
         ones(1, numFrames); ...       %  21
         data_reshaped(44:48, :); ...  %  22 -->  26
         zeros(6, numFrames);];
end