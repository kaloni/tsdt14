clear all
close all

% Parameters
num_samples = 2^10;
nInputs = 10000;
sampling_frequency = 10;
max_time = num_samples / sampling_frequency;
t = linspace(0,max_time,num_samples);

s = randn(num_samples, nInputs); % noise
cut_off = 0.1;
[h_b, h_a] = butter(7,2*cut_off); %sinc(t / sampling_frequency);
x = filter(h_b,h_a,s); % filtered noise

%% functions
omega = 1/2; % carrier frequency
%f1 = @(x) x .* repmat(((-1).^(0:(size(x,1)-1))), size(x,2));
%f2 = @(x) x.*(mod(0:(length(x)-1)) == 0);

%y1 = f1(x);
%y2 = f2(x);
y1 = x;
y2 = x;
for k = 1:size(x,2)
    y1(:,k) = x(:,k) .* ((-1).^(0:(size(x,1)-1)))';
end
for k = 1:size(x,2)
    y2(:,k) = x(:,k) .* (mod( 0:(size(x,1)-1), 2) == 0)';
end


disp('ok')

%% PSD in theory
theta = linspace(-1/2,1/2,num_samples);

fc = ceil((num_samples)*omega);

RX = zeros(num_samples,1);
RX(1:ceil(num_samples*cut_off)) = 1;
RX(end:-1:end-floor(num_samples*cut_off)) = 1;
rx = ifft(RX);
rx = abs(rx);

R1 = circshift(RX,[fc,0])/4 + circshift(RX,[-fc,0])/4;

R2 = RX/4 +  circshift(RX,[fc,0])/16 +  circshift(RX,[-fc,0])/16;

%figure(1);

%% PSD estimations
PX = psd(x);
P1 = psd(y1);
P2 = psd(y2);

theta = linspace(-0.5,0.5,num_samples);
figure(1);
title('estimated psds')
subplot(3,1,1); hold on; plot(theta,fftshift(PX),'c'); plot(theta,fftshift(RX),'r'); title('Input');
subplot(3,1,2); hold on; plot(theta,fftshift(P1),'c'); plot(theta,fftshift(R1),'r'); title('Alternating');
subplot(3,1,3); hold on; plot(theta,fftshift(P2),'c'); plot(theta,fftshift(R2),'r'); title('Decimation');
%hold on
%subplot(3,1,1); plot(theta,fftshift(RX),'r'); title('Input');
%hold on
%subplot(3,1,2); plot(theta,fftshift(R1),'r'); title('Alternating');
%hold on
%subplot(3,1,3); plot(theta,fftshift(R2),'r'); title('Decimation');
