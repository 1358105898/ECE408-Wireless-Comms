function [chan] = generateRayleighChannel(N, f_m)
%generateRayleighChannel generates a fading Rayleigh channel based on Fig.
%5.24 in Rappaport

%%
% (1) Specify the number of frequency domain points (N) used to represent
% sqrt(S_E(f)) and the maximum Doppler freqyency shift (f_m)
%	N   = number of frequency domain points, generally power of 2
%   f_m = maximum Doppler frequency shift

%%
% (2) Compute the frequency spacing between adjacent spectral lines. This
% defines the time duration of a fading waveform.
delta_f = 2*f_m/(N - 1);
freq_domain = (-f_m:delta_f:f_m).';

%%
% (3) Generate complex Gaussian random variables for each of the N/2
% positive frequency components of the noise source.
gaussian_randvars_pos = randn(ceil(N/2), 2) + 1j*randn(ceil(N/2), 2);

%%
% (4) Construct the negative frequency components of the noise source by
% conjugating positive frequency values and assigning these at negative
% frequency values.
switch mod(N, 2)
    case 0      % N is even
        gaussian_randvars = [flipud(conj(gaussian_randvars_pos)); ...
                             gaussian_randvars_pos];
    case 1      % N is odd
        gaussian_randvars = [flipud(conj(gaussian_randvars_pos(2:end, :))); ...
                             gaussian_randvars_pos];
end

%%
% (5) Multiply the in-phase and quadrature noise sources by the fading
% spectrum sqrt(S_E(f)).
%   Eq. 5.76 defines the Clarke Spectrum S_E_f; we'll be performing
%   baseband operations, so let f_c = 0. As well, to avoid going to inf at
%   f = f_m, we instead sample at f = 0.999999f_m
S_E_f = (1.5/(pi * f_m)) ./ sqrt(1 - ((freq_domain)/f_m).^2);
S_E_f([1, end]) = (1.5/(pi * f_m)) ./ sqrt(1 - ((freq_domain([1, end])*.999999)/f_m).^2);

ifft_in = gaussian_randvars .* sqrt(S_E_f);

%%
% (6) Perform an IFFT on the resulting frequency domain signals from the
% in-phase and quadrature arms to get two N-length time series, and add the
% squares of each signal point in time to create an N-point time series.
ifft_out = ifft(ifftshift(ifft_in), N);
sum_terms = sum(ifft_out.^2, 2);

%%
% (7) Take the square root of the sum to obtain an N-point time series of a
% simulated Rayleigh fading signal with the proper Doppler spread and time
% correlation.
chan = sqrt(sum_terms);
end

