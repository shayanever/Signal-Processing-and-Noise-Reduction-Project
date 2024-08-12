# Signal-Processing-and-Noise-Reduction-Project
A MATLAB-based project exploring various noise reduction techniques for audio signal processing. It includes scripts for applying and comparing filters like FIR, IIR, LMS, and Wiener to improve Signal-to-Noise Ratio (SNR) across different noise types, with visualizations and detailed analysis.


## Overview

This repository contains MATLAB scripts for a signal processing project focused on enhancing noisy audio signals. The project utilizes various filtering techniques to improve the Signal-to-Noise Ratio (SNR) and offers comprehensive visualizations to analyze the effectiveness of each method.

## Repository Structure

- **`app.m`**: 
  - A MATLAB script that applies multiple filters to a noisy signal, including Low-Pass, High-Pass, Band-Pass, Moving Average, and Wiener filters. It provides both time and frequency domain visualizations and features a simple GUI for audio playback.

- **`comparison.m`**: 
  - A script for loading a noisy audio signal and applying different filters. It generates visualizations of the signal in the time domain, frequency domain, and through spectrograms to compare the effects of each filter.

- **`snr_comp.m`**: 
  - This script evaluates the SNR improvements achieved by various noise reduction methods (FIR, IIR, LMS, RLS, Wiener, etc.) on a clean speech
