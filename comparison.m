close all;
clear all;
clc;
% Load the .wav file
%[noisy_signal, fs] = audioread('C:\Users\Shayan\Desktop\dsp project\New folder\three\create noisy audio\noisy_signal_complex.wav'); 
[noisy_signal, fs] = audioread('C:\Users\Shayan\Desktop\dsp project\New folder\three\noisy_datasetworking)_kagle\BWAVN\19-227-0046.wav'); 


% Generate a time vector based on the samplng frequenc
time = (0:length(noisy_signal)-1) / fs;

%% Apply Filters
% Low-Pass Filter
lpf_signal = lowPassFilter(noisy_signal, fs);

% High-Pass Filter
hpf_signal = highPassFilter(noisy_signal, fs);

% Band-Pass Filter
bpf_signal = bandPassFilter(noisy_signal, fs);

% Moving Average Filter
ma_signal = movingAverageFilter(noisy_signal, fs);

% Wiener Filter
wiener_signal = wienerFilter1D(noisy_signal, var(noisy_signal) * 0.01);

%% Time Domain Plots
figure;
subplot(6, 1, 1);
plot(time, noisy_signal);
title('Original Noisy Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6, 1, 2);
plot(time, lpf_signal);
title('Low-Pass Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6, 1, 3);
plot(time, hpf_signal);
title('High-Pass Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6, 1, 4);
plot(time, bpf_signal);
title('Band-Pass Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6, 1, 5);
plot(time, ma_signal);
title('Moving Average Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6, 1, 6);
plot(time, wiener_signal);
title('Wiener Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

%% Frequency Domain Plot
figure;
subplot(6, 1, 1);
fft_noisy = abs(fftshift(fft(noisy_signal)));
f = linspace(-fs/2, fs/2, length(fft_noisy));
plot(f, fft_noisy);
title('Original Noisy Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 2);
fft_lpf = abs(fftshift(fft(lpf_signal)));
plot(f, fft_lpf);
title('Low-Pass Filtered Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 3);
fft_hpf = abs(fftshift(fft(hpf_signal)));
plot(f, fft_hpf);
title('High-Pass Filtered Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 4);
fft_bpf = abs(fftshift(fft(bpf_signal)));
plot(f, fft_bpf);
title('Band-Pass Filtered Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 5);
fft_ma = abs(fftshift(fft(ma_signal)));
plot(f, fft_ma);
title('Moving Average Filtered Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 6);
fft_wiener = abs(fftshift(fft(wiener_signal)));
plot(f, fft_wiener);
title('Wiener Filtered Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% Spectrogram Plot
figure;
subplot(6, 1, 1);
spectrogram(noisy_signal, 256, [], [], fs, 'yaxis');
title('Original Noisy Signal - Spectrogram');

subplot(6, 1, 2);
spectrogram(lpf_signal, 256, [], [], fs, 'yaxis');
title('Low-Pass Filtered Signal - Spectrogram');

subplot(6, 1, 3);
spectrogram(hpf_signal, 256, [], [], fs, 'yaxis');
title('High-Pass Filtered Signal - Spectrogram');

subplot(6, 1, 4);
spectrogram(bpf_signal, 256, [], [], fs, 'yaxis');
title('Band-Pass Filtered Signal - Spectrogram');

subplot(6, 1, 5);
spectrogram(ma_signal, 256, [], [], fs, 'yaxis');
title('Moving Average Filtered Signal - Spectrogram');

subplot(6, 1, 6);
spectrogram(wiener_signal, 256, [], [], fs, 'yaxis');
title('Wiener Filtered Signal - Spectrogram');

%% FFT Magnitude Plots
figure;
subplot(6, 1, 1);
plot(f, fft_noisy);
title('Original Noisy Signal - FFT Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 2);
plot(f, fft_lpf);
title('Low-Pass Filtered Signal - FFT Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 3);
plot(f, fft_hpf);
title('High-Pass Filtered Signal - FFT Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 4);
plot(f, fft_bpf);
title('Band-Pass Filtered Signal - FFT Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 5);
plot(f, fft_ma);
title('Moving Average Filtered Signal - FFT Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(6, 1, 6);
plot(f, fft_wiener);
title('Wiener Filtered Signal - FFT Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% GUI for Audio Playback
% Create a figure for the GUI
hFig = figure('Name', 'Audio Player', 'NumberTitle', 'off');

uicontrol('Style', 'pushbutton', 'String', 'Play Original', ...
    'Position', [20 180 120 30], 'Callback', @(~,~) playAudio(noisy_signal, fs));
uicontrol('Style', 'pushbutton', 'String', 'Play Low-Pass Filtered', ...
    'Position', [20 140 120 30], 'Callback', @(~,~) playAudio(lpf_signal, fs));
uicontrol('Style', 'pushbutton', 'String', 'Play High-Pass Filtered', ...
    'Position', [20 100 120 30], 'Callback', @(~,~) playAudio(hpf_signal, fs));
uicontrol('Style', 'pushbutton', 'String', 'Play Band-Pass Filtered', ...
    'Position', [20 60 120 30], 'Callback', @(~,~) playAudio(bpf_signal, fs));
uicontrol('Style', 'pushbutton', 'String', 'Play Moving Average Filtered', ...
    'Position', [20 20 120 30], 'Callback', @(~,~) playAudio(ma_signal, fs));
uicontrol('Style', 'pushbutton', 'String', 'Play Wiener Filtered', ...
    'Position', [150 100 120 30], 'Callback', @(~,~) playAudio(wiener_signal, fs));

uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
    'Position', [150 60 80 30], 'Callback', @(~,~) pauseAudio());


uicontrol('Style', 'pushbutton', 'String', 'Resume', ...
    'Position', [150 20 80 30], 'Callback', @(~,~) resumeAudio());

%% Function for  Wiener Filter
function filtered_signal = wienerFilter1D(signal, noise_variance)
    signal_length = length(signal);
    
    % Estimate the signals power spectral density
    [pxx, f] = pwelch(signal, [], [], [], 'twosided');
    
    % Estimate the noise power spectral density 
    noise_psd = noise_variance * ones(size(pxx));
    
    % Calculate the Wiener filter in the frequency domain
    H = pxx ./ (pxx + noise_psd);
    
    % Interpolate the filter to match the length of the FFT of the signa
    H_interp = interp1(linspace(0, 1, length(H)), H, linspace(0, 1, signal_length));
    
    % Apply the Wiener filter in the frequency domain
    signal_fft = fft(signal);
    filtered_signal_fft = signal_fft .* H_interp';
    
    % Convert back to the time domain
    filtered_signal = real(ifft(filtered_signal_fft));
end

%% Additional Filter Function
function filtered_signal = lowPassFilter(signal, fs)
    cutoff_freq = 0.1; 
    order = 5; 
    [b, a] = butter(order, cutoff_freq, 'low');
    filtered_signal = filtfilt(b, a, signal);
end

function filtered_signal = highPassFilter(signal, fs)
    cutoff_freq = 0.3;
    order = 5; 
    [b, a] = butter(order, cutoff_freq, 'high');
    filtered_signal = filtfilt(b, a, signal);
end

function filtered_signal = bandPassFilter(signal, fs)
    low_cutoff = 0.2; 
    high_cutoff = 0.4;
    order = 6; 
    [b, a] = butter(order, [low_cutoff high_cutoff], 'bandpass');
    filtered_signal = filtfilt(b, a, signal);
end

function filtered_signal = movingAverageFilter(signal, fs)
    window_size = round(0.001 * fs); 
    filtered_signal = movmean(signal, window_size);
end

%% Helper Functions for Audio Playback
function playAudio(signal, fs)
    global player;
    player = audioplayer(signal, fs);
    play(player);
end

function pauseAudio()
    global player;
    if ~isempty(player)
        pause(player);
    end
end

function resumeAudio()
    global player;
    if ~isempty(player)
        resume(player);
    end
end
