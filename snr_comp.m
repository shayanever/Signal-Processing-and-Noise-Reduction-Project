clc;
clear all;
close all;

% Define Paths
cleanSpeechPath = 'c1.wav'; 
noiseFiles = {'white_noise.wav', 'pink_noise.wav', ... 
              'n3.wav', 'n1.wav'};

noiseNames = {'White Noise', 'Pink Noise', 'Bubble Noise', 'Engine Noise'};

% Load Clean Speech
[originalCleanSignal, fs_clean] = audioread(cleanSpeechPath);
originalCleanSignal = originalCleanSignal / max(abs(originalCleanSignal)); % Normalize clean signal

% Ensure cleanSignal is a column vector
originalCleanSignal = originalCleanSignal(:);

% Define SNR Levels
snrLevels = [-5, 1, 5]; 

% Preallocate results storage
snrResults = zeros(length(noiseFiles), length(snrLevels), 9); 
resultsTable = {}; 

% Loop Over Noise Files and SNR Levels
for noiseIdx = 1:length(noiseFiles)
    [noise, fs_noise] = audioread(noiseFiles{noiseIdx});
    if fs_clean ~= fs_noise
        noise = resample(noise, fs_clean, fs_noise);
    end
    noise = noise / max(abs(noise)); % Normalize noise

    % Ensure noise is a column vector
    noise = noise(:);

    % Ensure the noise length matches the clean signal length
    if length(noise) > length(originalCleanSignal)
        noise = noise(1:length(originalCleanSignal)); % Truncate noise to clean signal length
    elseif length(noise) < length(originalCleanSignal)
        % If noise is shorter, repeat to match length
        noise = repmat(noise, ceil(length(originalCleanSignal) / length(noise)), 1);
        noise = noise(1:length(originalCleanSignal)); % Ensure exact match
    end

    for snrIdx = 1:length(snrLevels)
        % Copy the clean signal for each SNR level and noise type
        cleanSignal = originalCleanSignal;
        
        % Mix Noise with Clean Signal at Desired SNR
        snr = snrLevels(snrIdx);
        signalPower = var(cleanSignal);
        noisePower = var(noise);
        scalingFactor = sqrt(signalPower / (10^(snr/10) * noisePower));

        % Ensure that both cleanSignal and noise are column vectors
        cleanSignal = cleanSignal(:);
        noise = noise(:);

        % Add noise to the clean 
        noisySignal = cleanSignal + scalingFactor * noise;
        
        % Ensure Length Match
        noisySignal = noisySignal(1:length(cleanSignal));

        % FIR Filter 
        orderFIR = 100; 
        cutoffFreqFIR = 3500; 

        bFIR = fir1(orderFIR, cutoffFreqFIR / (fs_clean / 2), 'low');
        filteredSignalFIR = filter(bFIR, 1, noisySignal);

        % IIR Filter Parameters
        orderIIR = 6; 
        cutoffFreqIIR = 3500; 

        [bIIR, aIIR] = butter(orderIIR, cutoffFreqIIR / (fs_clean / 2), 'low');
        filteredSignalIIR = filter(bIIR, aIIR, noisySignal);

        % LMS Adaptive Filtering
        muLMS = 0.01; 
        numWeights = 32; 
        lmsFilter = dsp.LMSFilter('Length', numWeights, 'StepSize', muLMS);
        [filteredSignalLMS, ~] = lmsFilter(noisySignal, cleanSignal);

        % RLS Adaptive Filtering
        lambdaRLS = 0.99; 
        deltaRLS = 0.1;  

        rlsFilter = dsp.RLSFilter('Length', numWeights, 'ForgettingFactor', lambdaRLS, 'InitialInverseCovariance', deltaRLS);
        [filteredSignalRLS, ~] = rlsFilter(noisySignal, cleanSignal);

        % Wiener Filtering
        nfft = 1024;
        [Pxx_noisy, F_noisy] = pwelch(noisySignal, hamming(nfft), nfft/2, nfft, fs_clean);
        [Pxx_noise, F_noise] = pwelch(noise, hamming(nfft), nfft/2, nfft, fs_clean);
        [S, F_stft, T] = stft(noisySignal, fs_clean, 'Window', hamming(nfft), 'OverlapLength', nfft/2, 'FFTLength', nfft);
        H_interp = interp1(F_noisy, max((Pxx_noisy - Pxx_noise) ./ Pxx_noisy, 0), F_stft, 'linear', 'extrap');
        H_interp = max(H_interp, 0); % Ensure no negative values
        S_filtered = S .* H_interp;
        filteredSignalWiener = istft(S_filtered, fs_clean, 'Window', hamming(nfft), 'OverlapLength', nfft/2, 'FFTLength', nfft);
        
        % Low-Pass Filter (LPF)
        cutoffFreqLPF = 3500; % Hz
        [bLPF, aLPF] = butter(4, cutoffFreqLPF / (fs_clean / 2), 'low');
        filteredSignalLPF = filter(bLPF, aLPF, noisySignal);

        % Moving Average Filter
        windowSize = 5; 
        filteredSignalMA = filter(ones(1, windowSize) / windowSize, 1, noisySignal);

        % Spectral Subtraction
        % Estimate noise spectrum from a known segment
        noiseEstimate = noise(1:nfft);
        [Pxx_noise_est, F_noise_est] = pwelch(noiseEstimate, hamming(nfft), nfft/2, nfft, fs_clean);

        % Convert Pxx_noise_est to match the length of the STFT segments
        Pxx_noise_est = interp1(F_noise_est, Pxx_noise_est, F_stft, 'linear', 'extrap');
        Pxx_noise_est = max(Pxx_noise_est, 0); % Ensure no negative values

        [S_noisy, F_noisy_SS, T_SS] = stft(noisySignal, fs_clean, 'Window', hamming(nfft), 'OverlapLength', nfft/2, 'FFTLength', nfft);
        H_ss = max(abs(S_noisy) - Pxx_noise_est, 0); % Spectral subtraction
        S_ss_smoothed = medfilt2(abs(S_noisy) - Pxx_noise_est, [5 5]); % Smoothing
        S_ss_filtered = S_noisy .* (S_ss_smoothed ./ abs(S_noisy));
        filteredSignalSpectralSubtraction = istft(S_ss_filtered, fs_clean, 'Window', hamming(nfft), 'OverlapLength', nfft/2, 'FFTLength', nfft);

        % Time-Frequency Analysis - STFT
        [S_TFA, F_TFA, T_TFA] = stft(noisySignal, fs_clean, 'Window', hamming(nfft), 'OverlapLength', nfft/2, 'FFTLength', nfft);
      
        S_TFA_filtered = S_TFA;
        threshold_TFA = 0.1 * max(abs(S_TFA(:)));
        S_TFA_filtered(abs(S_TFA) < threshold_TFA) = 0;
        filteredSignalTFA = istft(S_TFA_filtered, fs_clean, 'Window', hamming(nfft), 'OverlapLength', nfft/2, 'FFTLength', nfft);

        % Truncate to Match Lengths
        minLen = min([length(cleanSignal), length(filteredSignalFIR), length(filteredSignalIIR), ...
                      length(filteredSignalLMS), length(filteredSignalRLS), length(filteredSignalWiener), length(filteredSignalLPF), length(filteredSignalMA), ...
                      length(filteredSignalSpectralSubtraction), length(filteredSignalTFA)]);
        cleanSignal = cleanSignal(1:minLen);
        noisySignal = noisySignal(1:minLen);
        filteredSignalFIR = filteredSignalFIR(1:minLen);
        filteredSignalIIR = filteredSignalIIR(1:minLen);
        filteredSignalLMS = filteredSignalLMS(1:minLen);
        filteredSignalRLS = filteredSignalRLS(1:minLen);
        filteredSignalWiener = filteredSignalWiener(1:minLen);
        filteredSignalLPF = filteredSignalLPF(1:minLen);
        filteredSignalMA = filteredSignalMA(1:minLen);
        filteredSignalSpectralSubtraction = filteredSignalSpectralSubtraction(1:minLen);
        filteredSignalTFA = filteredSignalTFA(1:minLen);
        
        % Normalize signals before saving to avoid clipping
        normalizedFilteredSignalFIR = filteredSignalFIR / max(abs(filteredSignalFIR));
        normalizedFilteredSignalIIR = filteredSignalIIR / max(abs(filteredSignalIIR));
        normalizedFilteredSignalLMS = filteredSignalLMS / max(abs(filteredSignalLMS));
        normalizedFilteredSignalRLS = filteredSignalRLS / max(abs(filteredSignalRLS));
        normalizedFilteredSignalWiener = filteredSignalWiener / max(abs(filteredSignalWiener));
        normalizedFilteredSignalLPF = filteredSignalLPF / max(abs(filteredSignalLPF));
        normalizedFilteredSignalMA = filteredSignalMA / max(abs(filteredSignalMA));
        normalizedFilteredSignalSpectralSubtraction = filteredSignalSpectralSubtraction / max(abs(filteredSignalSpectralSubtraction));
        normalizedFilteredSignalTFA = filteredSignalTFA / max(abs(filteredSignalTFA));

        % Calculate SNR Improvement
        originalSNR = 10 * log10(var(cleanSignal) / var(noisySignal - cleanSignal));
        SNR_FIR = 10 * log10(var(cleanSignal) / var(filteredSignalFIR - cleanSignal));
        SNR_IIR = 10 * log10(var(cleanSignal) / var(filteredSignalIIR - cleanSignal));
        SNR_LMS = 10 * log10(var(cleanSignal) / var(filteredSignalLMS - cleanSignal));
        SNR_RLS = 10 * log10(var(cleanSignal) / var(filteredSignalRLS - cleanSignal));
        SNR_Wiener = 10 * log10(var(cleanSignal) / var(filteredSignalWiener - cleanSignal));
        SNR_LPF = 10 * log10(var(cleanSignal) / var(filteredSignalLPF - cleanSignal));
        SNR_MA = 10 * log10(var(cleanSignal) / var(filteredSignalMA - cleanSignal));
        SNR_SpectralSubtraction = 10 * log10(var(cleanSignal) / var(filteredSignalSpectralSubtraction - cleanSignal));
        SNR_TFA = 10 * log10(var(cleanSignal) / var(filteredSignalTFA - cleanSignal));

        % Store SNR improvements for plotting
        snrResults(noiseIdx, snrIdx, :) = [SNR_FIR, SNR_IIR, SNR_LMS, SNR_RLS, SNR_Wiener, SNR_LPF, SNR_MA, SNR_SpectralSubtraction, SNR_TFA];

        % Store all results in the resultsTable
        filterNames = {'FIR', 'IIR', 'LMS', 'RLS', 'Wiener', 'LPF', 'Moving Average', 'Spectral Subtraction', 'Time-Frequency Analysis'};
        snrValues = [SNR_FIR, SNR_IIR, SNR_LMS, SNR_RLS, SNR_Wiener, SNR_LPF, SNR_MA, SNR_SpectralSubtraction, SNR_TFA];
        for i = 1:length(filterNames)
            resultsTable{end+1, 1} = noiseNames{noiseIdx};
            resultsTable{end, 2} = snr;
            resultsTable{end, 3} = filterNames{i};
            resultsTable{end, 4} = snrValues(i);
        end
    end
end


resultsTable = cell2table(resultsTable, 'VariableNames', {'NoiseType', 'SNR', 'FilterType', 'SNRValue'});
disp(resultsTable);

% Plotting SNR improvements
methods = {'FIR', 'IIR', 'LMS', 'RLS', 'Wiener', 'LPF', 'Moving Average', 'Spectral Subtraction', 'Time-Frequency Analysis'};
for noiseIdx = 1:length(noiseFiles)
    figure;
    plot(snrLevels, squeeze(snrResults(noiseIdx, :, :)), '-o');
    title(sprintf('SNR Improvement for %s', noiseNames{noiseIdx}));
    xlabel('SNR (dB)');
    ylabel('SNR Improvement (dB)');
    legend(methods, 'Location', 'Best');
    grid on;
end

% Plot Time-Frequency Representations for comparison
for noiseIdx = 1:length(noiseFiles)
    figure;
    subplot(3,1,1);
    spectrogram(originalCleanSignal, hamming(nfft), nfft/2, nfft, fs_clean, 'yaxis');
    title('Original Clean Signal');
    subplot(3,1,2);
    spectrogram(noisySignal, hamming(nfft), nfft/2, nfft, fs_clean, 'yaxis');
    title(sprintf('Noisy Signal (%s at %d dB SNR)', noiseNames{noiseIdx}, snrLevels(snrIdx)));
    subplot(3,1,3);
    spectrogram(filteredSignalLMS, hamming(nfft), nfft/2, nfft, fs_clean, 'yaxis');
    title('Best Filtered Signal (LMS)');
end
