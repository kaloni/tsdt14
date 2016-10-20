function [ R ] = psd( x )
% Estimates the acf of a signal x
% Uses barlett estimate for the periodogram

f = fft(x);
R = abs(f).^2 / size(x,1);
R = mean(R,2);

%r = acf(x);
%R = fft(r);

end

