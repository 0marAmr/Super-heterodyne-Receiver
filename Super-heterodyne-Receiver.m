%% clear section
clear
clc
close all

%% message signal

filename = ["Short_QuranPalestine.wav", "Short_BBCArabic2.wav"];              % names of audio files
for i = 1 : 2
    [y, Fs] = audioread(filename(i));                                         % read audio and getting sample frequency
    size(y)
    y = sum(y.');                                                             % converte from stereo to mono sound
    y = y';                                                                   % converte from stereo to mono sound

    max_audio_length = 739200;                                                % max audio length
    y(end + max_audio_length - length(y) , 1) = 0;                            % adding zeros to short audios
    % sound(y,Fs)                                           % -+lay audio

    subplot(3,2,i)
    plot(y, 'b');
    title("(" + filename(i) + ") " + "Signal in Time Domain")
    xlabel("Time (s)")
    ylabel("Volts")
    grid on 
    ylim([-3 3])

    %% fft of message signal 
    Y = fftshift(fft(y));
    f_Y = -length(Y)/2:1:length(Y)/2-1;
    subplot(3,2,4 - mod(i,2))
    plot(f_Y, abs(Y), 'r');
    title("Frequency Spectrum of Signal")
    xlabel("Frequency (Hz)")
    ylabel("Magnitude")
    grid on 

    %% carrier signal
    n = (i-1);
    delta_f = 50*10^3;               % 50 KHz
    f_o = 100000;                    % 100 KHz
    fn = f_o + n*delta_f;            % carrier frequency
    Fs_carrier = 10*fn;
    t = (0:1/Fs_carrier:1-1/Fs_carrier)';
    carrier_signal = cos(2*pi*fn*t);

    %% AM modulation
    % modulating_signal = interp(y, Fs_carrier/length(y));
    [NUM, DEN] = numden(sym(length(carrier_signal)/length(y)));
    modulating_signal = resample(y, double(NUM), double(DEN));
    modulated_signal = carrier_signal.*modulating_signal; 
    MODULATED_SIGNAL = fftshift(fft(modulated_signal));

    f_MODULATED_SIGNAL = -length(MODULATED_SIGNAL)/2:1:length(MODULATED_SIGNAL)/2-1;
    subplot(3,2,6 - mod(i,2))
    plot(f_MODULATED_SIGNAL,abs(MODULATED_SIGNAL), 'g');
    title("Frequency Spectrum of Modulated Signal")
    xlabel("Frequency (Hz)")
    ylabel("Magnitude")
    grid on 
    xlim([-0.5*10^6 0.5*10^6])
end