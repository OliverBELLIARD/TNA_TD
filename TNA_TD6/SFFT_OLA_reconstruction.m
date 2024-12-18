close all; clear all; clc;

% Parameters
x = sin(2*pi*0.01*(0:999))'; % Example signal: a sine wave
M = 256; % Window size
% R = M/2; % 50% Overlap
R = M/4; % 75% Overlap
w = hamming(M); % Analysis window

% STFT Analysis
num_frames = floor((length(x) - M) / R) + 1;
X = zeros(M, num_frames);
for k = 0:num_frames-1
    x_seg = x(k*R + (1:M)); % Extract frame
    x_seg_windowed = x_seg .* w; % Apply window
    X(:,k+1) = fft(x_seg_windowed); % FFT
end

%% ISTFT and Overlap-Add Reconstruction
x_reconstructed = zeros(length(x), 1); % Initialize
w_synthesis = w; % Synthesis window (same as analysis)
for k = 0:num_frames-1
    x_ifft = ifft(X(:,k+1)); % Inverse FFT
    x_overlap = real(x_ifft) .* w_synthesis; % Apply window
    x_reconstructed(k*R + (1:M)) = x_reconstructed(k*R + (1:M)) + x_overlap;
end

% Plot results
figure;
subplot(2,1,1);
plot(x);
title('Original Signal');
xlabel('Samples');
ylabel('Amplitude');

subplot(2,1,2);
plot(x_reconstructed);
title('Reconstructed Signal using STFT and OLA');
xlabel('Samples');
ylabel('Amplitude');

%% Constant-overlap-add (COLA) constraint deconstruction
% Parameters
fs = 1000;               % Sampling frequency
t = 0:1/fs:5-1/fs;       % Time vector for 5 seconds
x = chirp(t,0,5,250);    % Test signal: Chirp signal from 0 Hz to 250 Hz

M = 512;                 % Window size
overlap = 50;            % Overlap percentage (50% in this case)
R = M * (1 - overlap/100); % Hop size

w = hamming(M);          % Window: Hamming
num_frames = floor((length(x) - M) / R) + 1;

%% Plotting setup
figure;

% Plot 1: Original Signal Sectioned into Buffers
subplot(5,1,1);
plot(t, x);
title('Complete Data Stream Sectioned into Buffers');
ylabel('Amplitude');
xlim([0, 5]);

% Plot 2: Windowed Buffers with Overlap
subplot(5,1,2);
buffered = zeros(M, num_frames);
for k = 0:num_frames-1
    idx = k*R + (1:M);
    buffered(:,k+1) = x(idx)' .* w; % Apply window
    hold on;
    plot(t(idx), w, 'Color', [0.5, 0.5, 0.5]); % Show windows
end
title('Windowing of Buffers with 50% Overlap');
ylabel('Amplitude');
ylim([0 1]);
xlim([0, 5]);

% Step 1: FFT, Process, and IFFT
processed_buffers = zeros(M, num_frames);
for k = 0:num_frames-1
    Xk = fft(buffered(:,k+1));   % Step 1: FFT
    Yk = Xk;                     % Step 2: Example signal manipulation (Identity)
    processed_buffers(:,k+1) = real(ifft(Yk)); % Step 3: IFFT
end

% Plot 3: Processed Buffers in Frequency Domain
subplot(5,1,3);
for k = 1:num_frames
    hold on;
    plot(t((k-1)*R + (1:M)), processed_buffers(:,k)); % Processed segments
end
title('Individual Processing of Each Windowed Signal Part');
ylabel('Amplitude');
xlim([0, 5]);

% Step 4: Overlap-Add Reconstruction
x_reconstructed = zeros(length(x), 1);
for k = 0:num_frames-1
    idx = k*R + (1:M);
    x_reconstructed(idx) = x_reconstructed(idx) + processed_buffers(:,k+1);
end

% Plot 4: Reconstructed Signal with Overlap-Add
subplot(5,1,4);
plot(t, x_reconstructed, 'b');
title('Overlap Add for Reconstruction into Stream Analysis');
ylabel('Amplitude');
xlim([0, 5]);

% Plot 5: Original vs Reconstructed Signal
subplot(5,1,5);
plot(t, x, 'k', t, x_reconstructed, 'r--');
legend('Original Signal', 'Reconstructed Signal');
title('Original vs Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 5]);
