clear; clc; close all;

% Parameters
N = 32; % Filter order
Nl = N - 1; % Filter length
wp = 0.43 * pi; % Passband edge frequency

% Design the 2-channel orthogonal filter bank
[h0, h1, g0, g1] = firpr2chfb(Nl, wp/pi);

% Tree-structured decomposition for analysis filters
% Analysis filters for an 8-channel bank
H = cell(3, 8); % Store filter coefficients at each level
H{1, 1} = h0; % First level
H{1, 2} = h1;

% Create filters for deeper levels
for level = 2:3 % Level 2 and 3
    for k = 1:2^(level-1)
        H{level, 2*k-1} = conv(H{level-1, k}, h0); % Lowpass branch
        H{level, 2*k} = conv(H{level-1, k}, h1);  % Highpass branch
    end
end

% Final set of 8 analysis filters
H8 = H(3, :); % Level 3 has 8 filters

% Generate the rectangular signal
signal = [ones(1, 200), zeros(1, 1000)];

% Decompose the signal using the 8-channel filter bank
channel_outputs = cell(1, 8);
for i = 1:8
    channel_outputs{i} = conv(signal, H8{i}); % Convolve signal with each filter
end

% Plot time responses of all channels
figure;
for i = 1:8
    subplot(8, 1, i);
    plot(channel_outputs{i});
    title(['Channel ', num2str(i), ' Time Response']);
    xlabel('Samples');
    ylabel('Amplitude');
end

% Check SBC Impulse Response
% Compute the overall synthesis response
% Create synthesis filters for 8-channel tree structure (similar to analysis)
G = cell(3, 8); % Store synthesis filters
G{1, 1} = g0; % First level
G{1, 2} = g1;

for level = 2:3
    for k = 1:2^(level-1)
        G{level, 2*k-1} = conv(G{level-1, k}, g0); % Lowpass branch
        G{level, 2*k} = conv(G{level-1, k}, g1);  % Highpass branch
    end
end

% Final synthesis filters
G8 = G(3, :);

% Compute the overall SBC impulse response
sbc_response = zeros(1, length(signal) + length(G8{1}) - 1);

for i = 1:8
    sbc_response = sbc_response + conv(channel_outputs{i}, G8{i});
end

% Plot SBC impulse response
figure;
stem(sbc_response, 'filled');
title('SBC Impulse Response');
xlabel('Samples');
ylabel('Amplitude');
grid on;
