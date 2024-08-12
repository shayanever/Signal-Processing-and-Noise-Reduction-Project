function noise_filtering_gui_v6
    
    initGUI();

    function initGUI()
        close all;

        % Create the GUI figur
        f = figure('Name', 'Noise Filtering GUI V6', 'NumberTitle', 'off', ...
                   'Position', [100, 100, 1400, 900]);

        % Create axe for time domain, frequency domain, and spectrogram plot
        ax1 = axes('Parent', f, 'Units', 'Normalized', 'Position', [0.1, 0.7, 0.35, 0.15]);
        ax2 = axes('Parent', f, 'Units', 'Normalized', 'Position', [0.55, 0.7, 0.35, 0.15]);
        ax3 = axes('Parent', f, 'Units', 'Normalized', 'Position', [0.1, 0.5, 0.35, 0.15]);
        ax4 = axes('Parent', f, 'Units', 'Normalized', 'Position', [0.55, 0.5, 0.35, 0.15]);
        ax5 = axes('Parent', f, 'Units', 'Normalized', 'Position', [0.1, 0.3, 0.35, 0.15]);
        ax6 = axes('Parent', f, 'Units', 'Normalized', 'Position', [0.55, 0.3, 0.35, 0.15]);

        % Create a button to load an audio file
        uicontrol('Style', 'pushbutton', 'String', 'Load Audio File', ...
                  'Units', 'Normalized', 'Position', [0.1, 0.1, 0.2, 0.08], ...
                  'Callback', @loadAudioFile);

        % Create button to play, pause, and resume original and filtered audio
        uicontrol('Style', 'pushbutton', 'String', 'Play Original Audio', ...
                  'Units', 'Normalized', 'Position', [0.35, 0.1, 0.2, 0.08], ...
                  'Callback', @playOriginalAudio);

        uicontrol('Style', 'pushbutton', 'String', 'Play Filtered Audio', ...
                  'Units', 'Normalized', 'Position', [0.6, 0.1, 0.2, 0.08], ...
                  'Callback', @playFilteredAudio);

        uicontrol('Style', 'pushbutton', 'String', 'Pause Audio', ...
                  'Units', 'Normalized', 'Position', [0.85, 0.1, 0.2, 0.08], ...
                  'Callback', @pauseAudio);

        % Create text boxes to dsplay detected noise type, applied filters, and filter parameters
        noiseText = uicontrol('Style', 'text', 'String', 'Noise Type: N/A', ...
                              'Units', 'Normalized', 'Position', [0.1, 0.94, 0.8, 0.04], ...
                              'HorizontalAlignment', 'center');
        filterText = uicontrol('Style', 'text', 'String', 'Applied Filters: N/A', ...
                               'Units', 'Normalized', 'Position', [0.1, 0.9, 0.8, 0.04], ...
                               'HorizontalAlignment', 'center');
        filterParamsText = uicontrol('Style', 'text', 'String', 'Fs: N/A, Fstop: N/A, Fpass: N/A, Order: N/A', ...
                                     'Units', 'Normalized', 'Position', [0.1, 0.86, 0.8, 0.04], ...
                                     'HorizontalAlignment', 'center');
        metricsText = uicontrol('Style', 'text', 'String', 'Evaluation Metrics: N/A', ...
                                'Units', 'Normalized', 'Position', [0.1, 0.82, 0.8, 0.04], ...
                                'HorizontalAlignment', 'center');

        % Initialize variables for storing audio data
        original_audio = [];
        filtered_audio = [];
        original_player = [];
        filtered_player = [];
        Fs = 0;
        Fstop = NaN;
        Fpass = NaN;
        Order = NaN;

        % Callback function to load an audio file
        function loadAudioFile(~, ~)
            % Clear previous data
            original_audio = [];
            filtered_audio = [];
            original_player = [];
            filtered_player = [];
            Fs = 0;
            Fstop = NaN;
            Fpass = NaN;
            Order = NaN;

            % Reset GUI
            cla(ax1);
            cla(ax2);
            cla(ax3);
            cla(ax4);
            cla(ax5);
            cla(ax6);

            set(noiseText, 'String', 'Noise Type: N/A');
            set(filterText, 'String', 'Applied Filters: N/A');
            set(filterParamsText, 'String', 'Fs: N/A, Fstop: N/A, Fpass: N/A, Order: N/A');
            set(metricsText, 'String', 'Evaluation Metrics: N/A');

            [file, path] = uigetfile('*.wav', 'Select an Audio File');
            if isequal(file, 0)
                disp('No file selected');
                return;
            end

            % Load the selected audio file
            [original_audio, Fs] = audioread(fullfile(path, file));
            original_audio = original_audio(:,1); % Use one channel if stereo

            % Detect noise type and apply filtering 
            detected_noise_type = detectNoiseTypeEnhanced(original_audio, Fs);
            [filtered_audio, applied_filters, Fstop, Fpass, Order] = applyFilteringPipelineEnhanced(original_audio, detected_noise_type, Fs);

            % Evaluate the filtering performance
            evaluation_metrics = evaluateFiltering(original_audio, filtered_audio);

            % Create audioplayer objects for original and filtered audio
            original_player = audioplayer(original_audio, Fs);
            filtered_player = audioplayer(filtered_audio, Fs);

            % Update text boxes with detected noise type, applied filters, and filter parameters
            set(noiseText, 'String', ['Noise Type: ', detected_noise_type]);
            set(filterText, 'String', ['Applied Filters: ', applied_filters]);
            set(filterParamsText, 'String', sprintf('Fs: %d Hz, Fstop: %.2f Hz, Fpass: %.2f Hz, Order: %d', Fs, Fstop, Fpass, Order));
            set(metricsText, 'String', sprintf('SNR Improvement: %.2f dB, MSE: %.4f', evaluation_metrics.SNR_Improvement, evaluation_metrics.MSE));

            % Update plots sequentially
            updatePlots();
        end

        function updatePlots()
            % Clear and update time domain plots
            updateTimeDomainPlot(ax1, original_audio, Fs, 'Original Audio Signal (Time Domain)');
            updateTimeDomainPlot(ax2, filtered_audio, Fs, 'Filtered Audio Signal (Time Domain)');

            % Clear and update frequency domain plots
            updateFrequencyDomainPlot(ax3, original_audio, Fs, 'Original Audio Frequency Domain');
            updateFrequencyDomainPlot(ax4, filtered_audio, Fs, 'Filtered Audio Frequency Domain');

            % Clear and update spectrograms
            updateSpectrogramPlot(ax5, original_audio, Fs, 'Original Audio Spectrogram');
            updateSpectrogramPlot(ax6, filtered_audio, Fs, 'Filtered Audio Spectrogram');
        end

        function updateTimeDomainPlot(ax, audio, Fs, titleText)
            axes(ax);
            cla(ax);
            plot((0:length(audio)-1)/Fs, audio);
            title(titleText);
            xlabel('Time (s)');
            ylabel('Amplitude');
            xlim([0 length(audio)/Fs]);
            ylim([-1 1]);  
            drawnow;
        end

        function updateFrequencyDomainPlot(ax, audio, Fs, titleText)
            axes(ax);
            cla(ax);
            N = length(audio);
            Y = fft(audio);
            P2 = abs(Y/N);
            P1 = P2(1:floor(N/2)+1);  
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(floor(N/2)))/N; 
            plot(f, P1);
            xlabel('Frequency (Hz)');
            ylabel('Magnitude');
            title(titleText);
            xlim([0 max(f)]);
            if max(P1) > 0
                ylim([0 max(P1)*1.1]);  
            end
            drawnow;
        end

        function updateSpectrogramPlot(ax, audio, Fs, titleText)
            axes(ax);
            cla(ax);
            audio = checkAndHandleNaNInf(audio); % Ensure no NaN or Inf values
            spectrogram(audio, 256, [], [], Fs, 'yaxis');
            title(titleText);
            drawnow;
        end

        % Callback function to play original audio
        function playOriginalAudio(~, ~)
            if isempty(original_player)
                errordlg('No audio file loaded!', 'Error');
                return;
            end
            play(original_player);
        end

        % Calback function to play filtered audio
        function playFilteredAudio(~, ~)
            if isempty(filtered_player)
                errordlg('No audio file loaded or filtered!', 'Error');
                return;
            end
            play(filtered_player);
        end

        % Callback function to pause audio
        function pauseAudio(~, ~)
            if ~isempty(original_player) && isplaying(original_player)
                pause(original_player);
            end
            if ~isempty(filtered_player) && isplaying(filtered_player)
                pause(filtered_player);
            end
        end
    end
end

% Enhanced noise type detection with statistical analysis
function noise_type = detectNoiseTypeEnhanced(audio_signal, Fs)
    % Compute the power spectral density (PSD)
    [Pxx, F] = periodogram(audio_signal, [], [], Fs);

    % Analyze statistical properties
    mean_psd = mean(Pxx);
    variance_psd = var(Pxx);
    skewness_psd = skewness(Pxx);
    kurtosis_psd = kurtosis(Pxx);

    % Additional characteristics
    spectral_flatness = geomean(Pxx) / mean(Pxx); % Measure of flatnes
    max_psd_freq = F(Pxx == max(Pxx)); % Frequency with maximum powe

    % Thresholds for different noise types
    white_noise_threshold = 0.15;
    hum_frequency = 50; 
    gaussian_noise_threshold = 0.10;
    pink_noise_threshold = 0.05;
    spectral_flatness_threshold = 0.9;

    % Decision making based on multiple features
    if variance_psd > white_noise_threshold && spectral_flatness > spectral_flatness_threshold
        noise_type = 'White Noise';
    elseif any(abs(max_psd_freq - hum_frequency) < 5) 
        noise_type = 'Hum Noise';
    elseif variance_psd > gaussian_noise_threshold && kurtosis_psd < 3
        noise_type = 'Gaussian Noise';
    elseif mean(Pxx(F < 200)) > pink_noise_threshold
        noise_type = 'Pink Noise';
    else
        noise_type = 'Complex Noise';
    end
end

% Improved filtering pipeline with adaptive filters
function [filtered_audio, applied_filters, Fstop, Fpass, Order] = applyFilteringPipelineEnhanced(audio_signal, noise_type, Fs)
  
    filtered_audio = audio_signal;
    applied_filters = 'No filtering applied';
    Fstop = NaN;
    Fpass = NaN;
    Order = NaN;

    switch noise_type
        case 'White Noise'
            % Apply adaptive noise cancellation
            filtered_audio = applyAdaptiveFilter(audio_signal, Fs);
            applied_filters = 'Adaptive Filter (ANC), Spectral Subtraction, Wiener Filter';
            
            % Follow up with spectral subtraction and Wiener filtering
            filtered_audio = applySpectralSubtraction(filtered_audio, Fs);
            filtered_audio = applyWienerFilter(filtered_audio, Fs);

        case 'Hum Noise'
            % Notch filtering at 50 Hz or 60 Hz
            filtered_audio = applyNotchFilter(audio_signal, Fs, 50, 4);
            applied_filters = 'Notch Filter, Adaptive Filter (ANC), Wiener Filter';
            
            % Folow up with adaptive filtering and Wiener filtering
            filtered_audio = applyAdaptiveFilter(filtered_audio, Fs);
            filtered_audio = applyWienerFilter(filtered_audio, Fs);

        case 'Gaussian Noise'
            % Band-pass filtering and adaptive filtering
            filtered_audio = applyBandPassFilter(audio_signal, Fs, [500, 3000], [450, 3100], 6);
            applied_filters = 'Band-Pass Filter, Adaptive Filter (ANC), Wiener Filter';
            
            % Follow up with adaptive filtering and Wiener filtering
            filtered_audio = applyAdaptiveFilter(filtered_audio, Fs);
            filtered_audio = applyWienerFilter(filtered_audio, Fs);

        case 'Pink Noise'
            % Band-stop filtering followed by adaptive filtering
            filtered_audio = applyBandStopFilter(audio_signal, Fs, [50, 200], [25, 225], 6);
            applied_filters = 'Band-Stop Filter, Adaptive Filter (ANC), Wiener Filter';
            
            % Follow up with adaptive filtering and Wiener filtering
            filtered_audio = applyAdaptiveFilter(filtered_audio, Fs);
            filtered_audio = applyWienerFilter(filtered_audio, Fs);

        case 'Complex Noise'
            % Multi-stage filtering for complex noise
            filtered_audio = applyAdaptiveFilter(audio_signal, Fs);
            filtered_audio = applySpectralSubtraction(filtered_audio, Fs);
            filtered_audio = applyWienerFilter(filtered_audio, Fs);
            applied_filters = 'Adaptive Filter (ANC), Spectral Subtraction, Wiener Filter';

        otherwise
            % Default to no filtering
            filtered_audio = audio_signal;
            applied_filters = 'No filtering applied';
    end

    % Ensure the filtered audio is the same length as the original audio
    if length(filtered_audio) ~= length(audio_signal)
        warning('Filtered audio length does not match original audio length. Adjusting size.');
        filtered_audio = resizeSignal(filtered_audio, length(audio_signal));
    end
    
    % Ensure no NaN or Inf values remain in the signa
    filtered_audio = checkAndHandleNaNInf(filtered_audio);
    
    % Avoid eensure filtered audio is not just zero 
    if all(filtered_audio == 0)
        warning('Filtered audio is entirely zeros. Reverting to original signal.');
        filtered_audio = audio_signal;
    end
end

% Adaptive noise cancellation filter (LMS Filter)
function filtered_audio = applyAdaptiveFilter(audio_signal, Fs)
    % Create a reference noise signal (for simulation purposes)
    % In practice, this should be a known noise reference.
    noise_reference = randn(size(audio_signal)) * 0.1; % White noise reference
    
    % Adaptive filter (Least Mean Squares - LMS)
    lmsFilt = dsp.LMSFilter('Length', 32, 'StepSize', 0.01);
    
    % Apply the adaptive filter
    [~, filtered_audio] = lmsFilt(noise_reference, audio_signal);
end

% Function to handle NaN/Inf values in the audio signal
function audio_signal = checkAndHandleNaNInf(audio_signal)
    if any(isnan(audio_signal)) || any(isinf(audio_signal))
        % Replace NaN/Inf values with zeros
        audio_signal(isnan(audio_signal) | isinf(audio_signal)) = 0;
    end
end

% Existing filtering functions with added validation
function filtered_audio = applyLowPassFilter(audio_signal, Fs, Fpass, Fstop, Order)
    if Fpass <= 0 || Fpass >= Fs/2
        warning('Invalid low-pass filter frequencies. Skipping low-pass filtering.');
        filtered_audio = audio_signal;
        return;
    end
    
    % Normalize the frequencies
    Wn = Fpass / (Fs/2);
    
    % Design the low-pass filter
    [b_lp, a_lp] = butter(Order, Wn, 'low');
    
    % Apply the filter
    filtered_audio = filter(b_lp, a_lp, audio_signal);
end

function filtered_audio = applyNotchFilter(audio_signal, Fs, Fpass, Order)
    if Fpass <= 0 || Fpass >= Fs/2
        warning('Invalid notch filter frequency. Skipping notch filtering.');
        filtered_audio = audio_signal;
        return;
    end
    
    % Normalize the frequency
    Wn = Fpass / (Fs/2);
    
    % Design the notch filter
    [b_notch, a_notch] = iirnotch(Wn, Wn / 35);
    
    % Apply the filter
    filtered_audio = filter(b_notch, a_notch, audio_signal);
end

function filtered_audio = applyBandPassFilter(audio_signal, Fs, Fpass, Fstop, Order)
    if any(Fpass <= 0) || any(Fpass >= Fs/2)
        warning('Invalid bandpass filter frequencies. Skipping bandpass filtering.');
        filtered_audio = audio_signal;
        return;
    end
    
    % Normalize the frequencies
    Wn = Fpass / (Fs/2);
    
    % Design the bandpass filter
    [b_bp, a_bp] = butter(Order, Wn, 'bandpass');
    
    % Apply the filter
    filtered_audio = filter(b_bp, a_bp, audio_signal);
end

function filtered_audio = applyBandStopFilter(audio_signal, Fs, Fpass, Fstop, Order)
    if any(Fpass <= 0) || any(Fpass >= Fs/2)
        warning('Invalid band-stop filter frequencies. Skipping band-stop filtering.');
        filtered_audio = audio_signal;
        return;
    end
    
    % Normalize the frequencie
    Wn = Fpass / (Fs/2);
    
    % Design the band-stop filte
    [b_bs, a_bs] = butter(Order, Wn, 'stop');
    
    % Apply the filte
    filtered_audio = filter(b_bs, a_bs, audio_signal);
end

function filtered_audio = applyWienerFilter(audio_signal, Fs)
    noise_power = var(audio_signal) / 10; 
    signal_power = abs(fft(audio_signal)).^2;
    H_wiener = signal_power ./ (signal_power + noise_power);
    filtered_audio = real(ifft(H_wiener .* fft(audio_signal)));
end

% Updated function to handle NaN/Inf and perform spectral subtractio
function enhanced_audio = applySpectralSubtraction(filtered_audio, Fs)
    window_length = 512; 
    overlap_length = 256;
    nfft = 1024;

    % Compute STFT
    [S, F, T] = stft(filtered_audio, Fs, 'Window', hamming(window_length), 'OverlapLength', overlap_length, 'FFTLength', nfft);
    
    % Estimate noise spectrum using first few frame
    noise_estimate = mean(abs(S(:, 1:10)), 2);
    
    % Apply spectral subtractio
    S_clean = max(abs(S) - noise_estimate, 0); % Ensure non-negative spectrum
    S_clean = S_clean .* exp(1i * angle(S)); % Reconstruct complex spectrum
    
    % Invers STFT to get enhanced audio
    enhanced_audio = istft(S_clean, Fs, 'Window', hamming(window_length), 'OverlapLength', overlap_length, 'FFTLength', nfft);
    
    % Ensure the result is rea
    enhanced_audio = real(enhanced_audio);
end

% Function to evaluate filtering performanc
function metrics = evaluateFiltering(original_audio, filtered_audio)
    % Ensure size compatibility before evaluatio
    if length(original_audio) ~= length(filtered_audio)
        filtered_audio = resizeSignal(filtered_audio, length(original_audio));
    end
    
    % Calculate noise by subtracting filtered from original
    noise = original_audio - filtered_audio;
    
    % Avoid division by zero by checking for non-zero power in original signal
    if sum(original_audio.^2) > 0 && sum(noise.^2) > 0
        SNR_before = 10 * log10(sum(original_audio.^2) / sum(noise.^2));
        SNR_after = 10 * log10(sum(filtered_audio.^2) / sum(noise.^2));
        metrics.SNR_Improvement = SNR_after - SNR_before;
    else
        metrics.SNR_Improvement = NaN; 
    end

    % Compute Mean Squared Error (MSE)
    if sum((original_audio - filtered_audio).^2) > 0
        metrics.MSE = mean((original_audio - filtered_audio).^2);
    else
        metrics.MSE = 0; % If signals are identical
    end
end

% Function to resize signal to match the original length
function resized_signal = resizeSignal(signal, target_length)
    current_length = length(signal);
    if current_length > target_length
        resized_signal = signal(1:target_length);
    elseif current_length < target_length
        resized_signal = [signal; zeros(target_length - current_length, 1)];
    else
        resized_signal = signal;
    end
end
