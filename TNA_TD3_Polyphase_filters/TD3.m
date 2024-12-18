close all
clear
clc

%% 1. Design Filter for SAR around 100 dB for Resampling Factor 3
Fs = 1000;          % Sampling Frequency (Hz)
Fpass = 200;        % Passband Frequency (Hz)
Fstop = 300;        % Stopband Frequency (Hz)
Apass = 1;          % Passband Ripple (dB)
Astop = 100;        % Stopband Attenuation (dB)

% Normalize the frequencies
Wp = Fpass / (Fs / 2);  % Passband frequency normalized to Nyquist
Ws = Fstop / (Fs / 2);  % Stopband frequency normalized to Nyquist

% Use firpm to design the filter
n = 50;               % Filter order
b = firpm(n, [0 Wp Ws 1], [1 1 0 0]); % FIR filter design

% Visualize the Frequency Response of the Filter
figure;
freqz(b, 1, 1024, Fs);
title('Frequency Response of the Designed Filter');

%% 2. Polyphase Filter for Decimation by 3

M = 3;               % Decimation factor
nPolyphase = M;      % Number of phases (equal to M)
N = length(b);       % Length of the filter

% Adjust the filter length to be a multiple of M (pad with zeros)
if mod(N, M) ~= 0
    padding = M - mod(N, M);  % Calculate padding needed
    b = [b, zeros(1, padding)]; % Pad filter coefficients
end

% Reshape the filter coefficients into M phases
bPoly = reshape(b, nPolyphase, []); 

% Create an example signal (for illustration purposes)
Fs = 1000;           % Sampling Frequency (Hz)
t = 0:1/Fs:1-1/Fs;   % Time vector
x = cos(2*pi*50*t) + cos(2*pi*150*t); % Example signal

% Filter the signal using the polyphase filter
xPolyphase = zeros(1, length(x));

for m = 1:nPolyphase
    % Apply filter to every m-th sample (phase m) of the signal
    filtered = filter(bPoly(m,:), 1, x(m:M:end));
    xPolyphase((m:M:end)) = filtered;
end

% Visualize the results
figure;
subplot(2,1,1);
stem(t, x);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
stem(t, xPolyphase);
title('Signal after Polyphase Decimation by 3');
xlabel('Time (s)');
ylabel('Amplitude');

%% 3. Polyphase Filter for Interpolation by 3

L = 3;               % Interpolation factor
nPolyphase = L;      % Number of phases (equal to L)

% Adjust the filter length to be a multiple of L (pad with zeros)
if mod(N, L) ~= 0
    padding = L - mod(N, L);  % Calculate padding needed
    b = [b, zeros(1, padding)]; % Pad filter coefficients
end

% Reshape the filter coefficients into L phases
bPoly = reshape(b, nPolyphase, []); 

% Create an example signal (for illustration purposes)
x = cos(2*pi*50*t);  % Example signal

% Upsample the signal by a factor of L
xUpsampled = upsample(x, L);

% Filter the upsampled signal using the polyphase filter
xPolyphaseInterp = zeros(1, length(xUpsampled));
for m = 0:nPolyphase-1
    xPolyphaseInterp(L*m+1:end) = filter(bPoly(m+1,:), 1, xUpsampled(L*m+1:end));
end

% Visualize the results
figure;
subplot(2,1,1);
stem(t, x);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
t_up = (0:length(xUpsampled)-1) / (Fs * L); % New time vector after upsampling
stem(t_up, xPolyphaseInterp);
title('Signal after Polyphase Interpolation by 3');
xlabel('Time (s)');
ylabel('Amplitude');

%% 4. Polyphase Filter for Interpolation by 3/2

L = 3;                % Interpolation factor (Upsample by 3)
M = 2;                % Decimation factor (Downsample by 2)
nPolyphaseInterp = L; % Number of phases for interpolation
nPolyphaseDecim = M;  % Number of phases for decimation

% Adjust the filter length to be a multiple of L (pad with zeros)
if mod(N, L) ~= 0
    padding = L - mod(N, L);  % Calculate padding needed
    b = [b, zeros(1, padding)]; % Pad filter coefficients
end

% Reshape the filter coefficients into L phases for interpolation
bPolyInterp = reshape(b, nPolyphaseInterp, []);

% Adjust the filter length to be a multiple of M (pad with zeros)
if mod(N, M) ~= 0
    padding = M - mod(N, M);  % Calculate padding needed
    b = [b, zeros(1, padding)]; % Pad filter coefficients
end

% Reshape the filter coefficients into M phases for decimation
bPolyDecim = reshape(b, nPolyphaseDecim, []);    

% Create an example signal (for illustration purposes)
x = cos(2*pi*50*t);  % Example signal

% Upsample the signal by 3
xUpsampled = upsample(x, L);

% Apply polyphase filtering for interpolation by 3
xInterpPolyphase = zeros(1, length(xUpsampled));
for m = 0:nPolyphaseInterp-1
    xInterpPolyphase(L*m+1:end) = filter(bPolyInterp(m+1,:), 1, xUpsampled(L*m+1:end));
end

% Decimate the signal by 2
xFinal = downsample(xInterpPolyphase, M);

% Visualize the results
figure;
subplot(2,1,1);
stem(t, x);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
t_final = (0:length(xFinal)-1) / (Fs * L / M); % Time vector after both interpolation and decimation
stem(t_final, xFinal);
title('Signal after Polyphase Interpolation by 3/2');
xlabel('Time (s)');
ylabel('Amplitude');
