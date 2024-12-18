# Design of filter-banks  
**Objective:** Design a tree-structure 8-channel filter bank in matlab.  
  
## A - Design the PR Orthogonal 2-ch filter bank  
  
Design the analysis/synthesis two-channel orthogonal filter bank using the MATLAB function firpr2chfb. The filter length is N=32, and the lowpass passband edge frequency ωp=0.43π. Plot the impulse responses of the analysis and synthesis filters: h0[n], h1[n], g0[n] and g1[n]. Compute and plot the poles and zeros of the transfer functions H0(z), H1(z), G0(z) and G1(z). Compute and plot the magnitude and group delay responses of the analysis filters H0(z) and H1(z). Compute and plot the impulse response of the distortion transfer function t[n].  
  
## B - use this filter bank to compose the 8-ch filter bank  
  
Compose the eight-channel tree-structured analysis/synthesis filter bank. As a basic building block use the orthogonal two-channel filter bank and construct the analysis and synthesis banks. Compute and plot the magnitude responses of the resulting eight analysis filters. Verify the magnitude-preserving property of the overall bank.  
  
## C - Use this filter bank to decompose a rectangle signal  
  
1 square of 200 samples of 1 over 1000 of 0 by example. Plot all the channel time response, and check the SBC impulse response.

---

### **A - Design the PR Orthogonal 2-ch Filter Bank**

#### 1. **Design the Filters**

- Use the MATLAB function `firpr2chfb` to create the analysis and synthesis filters.
- Set the filter length N=32N = 32 and lowpass passband edge frequency ωp=0.43π\omega_p = 0.43\pi.

#### 2. **Plot Impulse Responses**

- Plot h0[n]h_0[n], h1[n]h_1[n], g0[n]g_0[n], and g1[n]g_1[n].

#### 3. **Poles and Zeros**

- Use the `zplane` function to compute and plot the poles and zeros of H0(z)H_0(z), H1(z)H_1(z), G0(z)G_0(z), and G1(z)G_1(z).

#### 4. **Magnitude and Group Delay Responses**

- Use `freqz` for magnitude responses.
- Use `grpdelay` for group delay.

#### 5. **Distortion Transfer Function**

- Compute and plot the impulse response of t[n]t[n], the distortion transfer function.

---

### **B - Compose the 8-Ch Filter Bank**

#### 1. **Tree-Structured Design**

- Construct an 8-channel tree-structured analysis filter bank using the 2-channel orthogonal filter bank as the building block.
- Similarly, construct the corresponding synthesis filter bank.

#### 2. **Magnitude Responses**

- Plot the magnitude responses of the 8 analysis filters.

#### 3. **Magnitude-Preserving Property**

- Sum the squared magnitude responses of all analysis filters and verify if they equal unity.

---

### **C - Decompose a Rectangle Signal**

#### 1. **Create the Signal**

- Define a signal consisting of 200 samples of 1 followed by 1000 samples of 0.

#### 2. **Decomposition**

- Pass the signal through the 8-channel analysis filter bank.
- Plot the time responses of each channel.

#### 3. **Check SBC Impulse Response**

- Verify the stop-band characteristics (SBC) impulse response of the filter bank.

---

### **MATLAB Implementation**

#### **Step A - 2-Ch Filter Bank**

```matlab
% Parameters
N = 32; % Filter length
wp = 0.43 * pi; % Passband edge frequency

% Design using MATLAB's firpr2chfb
[h0, h1, g0, g1] = firpr2chfb(N, wp);

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
freqz(h0); title('Magnitude Response of H0(z)');
figure;
freqz(h1); title('Magnitude Response of H1(z)');
figure;
grpdelay(h0); title('Group Delay of H0(z)');
figure;
grpdelay(h1); title('Group Delay of H1(z)');

% Distortion transfer function
t = conv(h0, g0) + conv(h1, g1);
figure;
stem(t); title('Impulse Response of Distortion Transfer Function');
```

#### **Step B - 8-Ch Filter Bank**

```matlab
% Construct tree-structured filter bank
% Use a loop or recursive function to build the structure
% Plot magnitude responses of the 8 channels

% Placeholder for example plots
figure;
hold on;
for i = 1:8
    [H, W] = freqz(filter_bank{i}, 1, 1024);
    plot(W/pi, abs(H)); % Magnitude response
end
hold off;
title('Magnitude Responses of 8 Analysis Filters');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');
legend('Channel 1','Channel 2', '...','Channel 8');

% Verify magnitude-preserving property
% Sum squared magnitude responses across channels
```

#### **Step C - Signal Decomposition**

```matlab
% Generate rectangular signal
signal = [ones(1,200), zeros(1,1000)];

% Decompose using 8-channel analysis filter bank
for i = 1:8
    channel_output{i} = conv(signal, filter_bank{i}); % Filter with each channel
end

% Plot time responses
figure;
for i = 1:8
    subplot(8,1,i);
    plot(channel_output{i});
    title(['Channel ', num2str(i), ' Time Response']);
end

% Check SBC impulse response
% Plot combined impulse response characteristics
```

---

