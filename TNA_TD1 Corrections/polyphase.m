clc
clearvars
close all

%% 1. Initialize Parameters
% These parameters define the input signal's properties, the resampling factors, and the filter settings.

FSin = 48000;       % Input signal's sampling frequency (Hz)
Up_Ratio = 3;       % Upsampling ratio (factor by which signal is upsampled)
Down_Ratio = 2;     % Downsampling ratio (factor by which signal is downsampled)
output_delay = [1 0];  % Output delay in polyphase fractional filtering (for fractional delay filter)
phase_length = 61;  % Length of the phase response in the FIR filter
coef_frac = [31 31]; % Fractional coefficient precision (not used directly in this script but likely related to the filter's precision)
nb_stages = length(Up_Ratio);  % Number of stages for resampling (used to determine the length of filter)

%% 2. Generate the Input Signal
% Create a cosine wave input signal with a frequency of 900 Hz.
t = 1/FSin:1/FSin:(max(Down_Ratio) * max(Up_Ratio) * 900) / FSin; % Time vector for 900 Hz signal, scaled by max up/down ratios
in = 2^24 * cos(2 * pi * 900 * t);  % Input signal (900 Hz cosine with high amplitude)

%% 3. Design the Filter for Polyphase Decimation and Fractional Interpolation
% Here, we design the filter coefficients for both decimation and interpolation. 

% Design the filter using the Parks-McClellan algorithm (firpm). The filter 
% length is calculated based on the upsampling and downsampling ratios.
h(1:(Up_Ratio * Down_Ratio * phase_length)) = firpm( ...
    (Up_Ratio * Down_Ratio * phase_length) - 1, ...   % Filter order
    [0, 0.9 / (max([Down_Ratio Up_Ratio])), 1 / (max([Down_Ratio Up_Ratio])), 1], ...   % Frequency bands
    [1, 1, 0, 0] ...   % Desired amplitude response: passband (1), stopband (0)
); 

% Plot the frequency response of the filter (to check the filter design)
freqz(h, 1, FSin / 20, FSin);

%% 4. Polyphase Decimation
% This step applies the polyphase decimation filter to the input signal.
% Polyphase decimation involves filtering the signal and then downsampling it by the decimation ratio.

out_decimated = polyphase_decimation(in, h, 3);  % Decimate the input signal by 3 using the polyphase filter

%% 5. Polyphase Fractional Interpolation
% This section performs fractional interpolation by combining upsampling and fractional delay filtering.

[out_fract, out_fract_direct] = polyphase_fractional(in, h, Up_Ratio, Down_Ratio, output_delay);  % Perform fractional interpolation

%% 6. Plot the Results
% Here we visualize the results of the decimated and interpolated signals.

% Generate x-axis indices for plotting
x = 1:500;          % Plotting the first 500 samples of the original input signal
x_dec = 1:3:500;    % Decimated signal, downsampled by 3

% Figure 2: Plot the original signal and the decimated signal
% The first subplot shows the original signal, and the second one shows the decimated signal.
figure(2)
plot(x, in(x), '-+', x_dec - 22.5, out_decimated(1:length(x_dec)), '-o');
% '-+' is used for the input signal markers, and '-o' for the decimated signal

% Figure 3: Compare the fractional interpolation results (both direct and using polyphase interpolation)
% The first subplot shows the direct fractional interpolation result, and the second one compares it with the polyphase fractional interpolation.
figure(3)
plot(x, out_fract_direct(x), '-+', x, out_fract(x + 1), '-o');
% '-+' is used for the direct interpolation result, and '-o' for the polyphase result

