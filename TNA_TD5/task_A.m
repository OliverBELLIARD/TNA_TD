clear all; clc; close all

% Parameters
N = 32; % Filter order
Nl = N-1; % Filter length
wp = 0.43 * pi; % Passband edge frequency

% Design using MATLAB's firpr2chfb
[h0, h1, g0, g1] = firpr2chfb(Nl, wp/pi);

% Impulse responses
figure;
subplot(2,2,1); stem(h0); title('Impulse Response h0[n]');
subplot(2,2,2); stem(h1); title('Impulse Response h1[n]');
subplot(2,2,3); stem(g0); title('Impulse Response g0[n]');
subplot(2,2,4); stem(g1); title('Impulse Response g1[n]');

% Poles and zeros
figure;
subplot(2,2,1); zplane(h0); title('Poles and Zeros of H0(z)');
subplot(2,2,2); zplane(h1); title('Poles and Zeros of H1(z)');
subplot(2,2,3); zplane(g0); title('Poles and Zeros of G0(z)');
subplot(2,2,4); zplane(g1); title('Poles and Zeros of G1(z)');

% Magnitude and group delay responses
figure;
freqz(h0);
hold on
freqz(h1);
title('Magnitude Responses of H0(z) & H1(z)');
legend('H0(z)', 'H1(z)');

figure;
grpdelay(h0);
hold on
grpdelay(h1);
title('Group Delays of H0(z) & H1(z)');
legend('H0(z)', 'H1(z)');

% Distortion transfer function
t = conv(h0, g0) + conv(h1, g1);
figure;
stem(t); title('Impulse Response of Distortion Transfer Function');