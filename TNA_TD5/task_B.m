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

% Magnitude responses of 8 analysis filters
figure;
hold on;
for i = 1:8
    [H8_resp, W] = freqz(H8{i}, 1, 1024);
    plot(W/pi, abs(H8_resp), 'DisplayName', ['Channel ', num2str(i)]);
end
title('Magnitude Responses of 8 Analysis Filters');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');
legend show;
grid on;

% Verify the magnitude-preserving property
% Sum squared magnitude responses across all channels
M = 0;
for i = 1:8
    [H8_resp, ~] = freqz(H8{i}, 1, 1024);
    M = M + abs(H8_resp).^2;
end

figure;
plot(W/pi, M);
title('Sum of Squared Magnitudes of Analysis Filters');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Sum of Magnitudes');
grid on;
