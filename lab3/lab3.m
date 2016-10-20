clear all
close all

% Parameters
num_samples = 2^10;
nInputs = 100;
sampling_frequency = 10;
max_time = num_samples / sampling_frequency;
t = linspace(0,max_time,num_samples);

s = randn(num_samples, nInputs); % noise
cut_off = 0.1;

[h_b, h_a] = butter(7,2*cut_off); %sinc(t / sampling_frequency);
x = filter(h_b,h_a,s); % filtered noise

%% functions
omega = 1/4; % carrier frequency
f1 = @(x) x.^2;
f2 = @(x) x.*(x>0);
f3 = @(x) x.*cos(2 * pi * omega * repmat((0:(size(x,1)-1))', 1, size(x,2)));

y1 = f1(x);
y2 = f2(x);
y3 = f3(x);
disp('ok')

%% PSD in theory
theta = linspace(-1/2,1/2,num_samples);
RX = zeros(num_samples,1);
RX(1:ceil(num_samples*cut_off)) = 1;
RX(end:-1:end-floor(num_samples*cut_off)) = 1;

rx = ifft(RX);
rx = abs(rx);

r1 = 2*rx.^2 + rx(1)^2;
R1 = abs(fft(r1));

r2 = rx/4 + (sqrt(rx(1)^2 - rx.^2) + rx.*asin(rx/rx(1)))/(2*pi);
R2 = abs(fft(r2));

fc = ceil((num_samples)*omega);
R3 = zeros(num_samples,1);
R3 = R3 + circshift(RX,[fc,0]);
R3 = R3 + circshift(RX,[-fc,0]);
R3 = R3 / 4;

PX = psd(x);
P1 = psd(y1);
P2 = psd(y2);
P3 = psd(y3);

figure(1)
A = 0.1;
subplot(2,2,1); hold on; plot(theta, log(A + fftshift(RX)),'r'); plot(theta, log(A + fftshift(PX)),'c'); title('Input');
subplot(2,2,2); hold on; plot(theta, log(A + fftshift(R1)),'r'); plot(theta, log(A + fftshift(P1)),'c'); title('Squarer');
subplot(2,2,3); hold on; plot(theta, log(A + fftshift(R2)),'r'); plot(theta, log(A + fftshift(P2)),'c'); title('Rectifier');
subplot(2,2,4); hold on; plot(theta, log(A + fftshift(R3)),'r'); plot(theta, log(A + fftshift(P3)),'c'); title('AM-SC');
% figure(1);
% title('theoretical psds')
% subplot(2,2,1); plot(theta,fftshift(RX)); title('Input PSD');
% subplot(2,2,2); plot(theta,fftshift(R1)); title('Squarer');
% subplot(2,2,3); plot(theta,fftshift(R2)); title('Rectifier');
% subplot(2,2,4); plot(theta,fftshift(R3)); title('AM-SC');
% 
% RX(1) = [];
% R1(1) = [];
% R2(1) = [];
% R3(1) = [];
% theta = linspace(-0.5,0.5,size(R1,1));
% figure(2);
% title('theoretical psds without dc')
% subplot(2,2,1); plot(theta,fftshift(RX)); title('Input');
% subplot(2,2,2); plot(theta,fftshift(R1)); title('Squarer');
% subplot(2,2,3); plot(theta,fftshift(R2)); title('Rectifier');
% subplot(2,2,4); plot(theta,fftshift(R3)); title('AM-SC');
% 
% %% PSD estimations
% PX = psd(x);
% P1 = psd(y1);
% P2 = psd(y2);
% P3 = psd(y3);
% 
% theta = linspace(-0.5,0.5,num_samples);
% figure(3);
% title('estimated psds')
% subplot(2,2,1); plot(theta,fftshift(PX)); title('Input')
% subplot(2,2,2); plot(theta,fftshift(P1)); title('Squarer');
% subplot(2,2,3); plot(theta,fftshift(P2)); title('Rectifier');
% subplot(2,2,4); plot(theta,fftshift(P3)); title('AM-SC');
% 
% P1 = P1(2:end);
% P2 = P2(2:end);
% P3 = P3(2:end);
% PX = PX(2:end);
% 
% theta = linspace(-0.5,0.5,size(P1,1));
% figure(4);
% title('estimated psds without dc')
% subplot(2,2,1); plot(theta,fftshift(PX)); title('Input')
% subplot(2,2,2); plot(theta,fftshift(P1)); title('Squarer');
% subplot(2,2,3); plot(theta,fftshift(P2)); title('Rectifier');
% subplot(2,2,4); plot(theta,fftshift(P3)); title('AM-SC');

%% PDF estimations (ampltitude hist)
figure(5)
subplot(2,2,1); hist(x(:),100); title('Histogram of input')
subplot(2,2,2); hist(y1(:),100); title('Histogram of squarer')
subplot(2,2,3); hist(y2(:),100); title('Histogram of rectifier')
subplot(2,2,4); hist(y3(:),100); title('Histogram of AM-SC')
