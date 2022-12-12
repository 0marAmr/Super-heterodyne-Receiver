%% fft for sin(2*pi*f*t)
clear
clc
close all

%% message signal
fm = 100;         % 10 KHz
Fs = 10*fm;
t = 0:1/Fs:1-1/Fs;
y = cos(2*pi*fm*t);
subplot(3,1,1)
plot(t,y, 'b');
title("Signal On Time Domain")
xlabel("Time (s)")
ylabel("Volts")
grid on 
ylim([-2 2])
xlim([0 2*1/fm])

%% fft of message signal 
Y = fftshift(fft(y));
f = -Fs/2:1:Fs/2-1;
subplot(3,1,2)
plot(f, abs(Y), 'r');
title("Frequency Spectrum of Signal")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
grid on 

%% carrier signal
n = 0;
delta_f = 50*10^3;              % 50 KHz
f_o = 1000;                     % 100 KHz
fn = f_o + n*delta_f;           % carrier frequency
Fs_carrier = 10*fn;
t = 0:1/Fs_carrier:1-1/Fs_carrier;
carrier_signal = cos(2*pi*fn*t);

%% AM modulation
modulating_signal = interp(y, Fs_carrier/length(y));
% [NUM, DEN] = numden(sym(Fs_carrier/length(y)));
% modulating_signal = resample(y, double(NUM), double(DEN));

% sound(y , 44100)
modulated_signal = carrier_signal.*modulating_signal; 
MODULATED_SIGNAL = fftshift(fft(modulated_signal));

f = -length(MODULATED_SIGNAL)/2:1:length(MODULATED_SIGNAL)/2-1;
subplot(3,1,3)
plot(f,abs(MODULATED_SIGNAL), 'k');
title("Frequency Spectrum of Modulated Signal")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
grid on 
axis 'auto xy'

%% test cos with high frequency
% clear
% clc
% close all
% Fs = 10000;
% t = 0 : 1/Fs : 1-1/Fs;  % Fs points
% Fm = 1000;
% y = cos(2*pi*Fm*t);
% plot(t, y)
% ylim([-2 2])
% xlim([0 0.01])