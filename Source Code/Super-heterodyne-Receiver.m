%% clear section
clear
clc
close all

%% message signal
filename = ["Short_BBCArabic2.wav", "Short_FM9090.wav", "Short_QuranPalestine.wav"];                      % names of audio files
sum_signals = 0;
for i = 1 : 3
    [y, Fs] = audioread(filename(i));                                         % read audio and getting sample frequency
    % size(y)
    y = sum(y.');                                                             % converte from stereo to mono sound
    y = y';                                                                   % converte from stereo to mono sound
    
    max_audio_length = 740544;                                                % max audio length
    y(end + max_audio_length - length(y) , 1) = 0;                            % adding zeros to short audios
    % sound(y,Fs)                                                             % play audio

    subplot(3,3,i)
    plot(y, 'b');
    title("(" + filename(i) + ") " + "Signal in Time Domain")
    xlabel("Time (s)")
    ylabel("Volts")
    grid on 
    ylim([-3 3])

    %% fft of message signal 
    Y = fft(y, length(y));
    f_Y = (-length(Y)/2:1:length(Y)/2-1)';
    subplot(3,3,3 + i)
    plot(f_Y*(Fs/length(Y)), fftshift(abs(Y)), 'r')
    title("Frequency Spectrum of Signal")
    xlabel("Frequency (Hz)")
    ylabel("Magnitude")
    grid on 

    %% carrier signal
    n = (i-1);
    delta_f = 50000;                 % 50 KHz
    f_o = 100000;                    % 100 KHz
    fn = f_o + n*delta_f;            % carrier frequency
    Fs_carrier = 10*Fs;              % 10*fn;
    t = (0:1/Fs_carrier:1-1/Fs_carrier)';
    carrier_signal = cos(2*pi*fn*t);

    %% AM modulation
    modulating_signal = interp(y, 10);
    carrier_signal(end + length(modulating_signal) - length(carrier_signal), 1) = 0;
    modulated_signal = carrier_signal.*modulating_signal; 
    MODULATED_SIGNAL = fftshift(fft(modulated_signal));
  
    sum_signals = sum_signals + MODULATED_SIGNAL;

    f_MODULATED_SIGNAL = (-length(MODULATED_SIGNAL)/2:1:length(MODULATED_SIGNAL)/2-1)';
    subplot(3,3,6 + i)
    plot(f_MODULATED_SIGNAL*Fs_carrier/length(MODULATED_SIGNAL),abs(MODULATED_SIGNAL), 'g')
    title("Frequency Spectrum of Modulated Signals")
    xlabel("Frequency (Hz)")
    ylabel("Magnitude")
    grid on 
    
end

f_MODULATED_SIGNAL = (-length(sum_signals)/2:1:length(sum_signals)/2-1)';
figure
% subplot(2,1,1)
plot(f_MODULATED_SIGNAL*Fs_carrier/length(sum_signals),abs(sum_signals), 'g')
title("Transmitted Signals")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
grid on 

%% RF Band Pass Filter
FpassLower = 141440;
FpassUpper = 158672;
sum_signals_BPF = bandpass(abs(sum_signals), [FpassLower, FpassUpper], 10*Fs);
SUM_SIGNALS_BPF = fftshift(fft(sum_signals_BPF));
f_MODULATED_SIGNAL = (-length(SUM_SIGNALS_BPF)/2:1:length(SUM_SIGNALS_BPF)/2-1)';
figure
plot(f_MODULATED_SIGNAL*Fs_carrier/length(SUM_SIGNALS_BPF), abs(SUM_SIGNALS_BPF), 'r')
title("After BPF")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
grid on
