
---

## **A- Verifying Constant Overlap-Add (COLA) Property**

### The **COLA property** guarantees that the overlap-add of shifted windows sums to a constant value, allowing perfect signal reconstruction without loss.

For the different windows:

1. **Rectangular Window at 0% overlap (Hop Size R=MR = M)**
    
    - No overlap between windows. Check if constant addition holds.
2. **Bartlett Window at 50% overlap (R≈M/2R \approx M/2)**
    
3. **Hamming Window at 50% overlap**
    
4. **Rectangular Window at 50% overlap**
    
5. **Hamming Window at 75% overlap (R=M/4R = M/4)**
    
6. **Blackman Window at 2/3 overlap (Hop Size = 1/3 MM)**
    

---

### MATLAB Code for A- COLA Check

```matlab
% Parameters
M_list = [32, 33, 64]; % Example window sizes
windows = {'rectwin', 'bartlett', 'hamming', 'rectwin', 'hamming', 'blackman'};
overlap_ratios = [0, 0.5, 0.5, 0.5, 0.75, 2/3];

figure;

for k = 1:length(windows)
    M = M_list(1); % Example M for demonstration (or choose appropriate M per window)
    overlap = overlap_ratios(k);
    R = round(M * (1 - overlap)); % Hop size

    % Define Window
    switch windows{k}
        case 'rectwin'
            w = rectwin(M);
        case 'bartlett'
            w = bartlett(M);
        case 'hamming'
            w = hamming(M);
        case 'blackman'
            w = blackman(M, 'periodic');
    end

    % Shift and Sum for COLA test
    num_shifts = 5; % Number of shifted windows to test
    cola_sum = zeros(M + R * (num_shifts-1), 1);
    for shift = 0:num_shifts-1
        cola_sum(shift*R + (1:M)) = cola_sum(shift*R + (1:M)) + w;
    end

    % Plot the result
    subplot(3, 2, k);
    plot(cola_sum);
    title(sprintf('%s Window (%.0f%% Overlap)', windows{k}, overlap*100));
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Constant Overlap-Add (COLA) Test');
```

**How to Analyze the Output**:

- For each subplot, verify whether the summed signal is **flat** (constant value).
- A flat curve confirms the COLA property.

---

## **B- Implement Short-Time Fourier Transform (STFT) and Verify OLA**

### Steps:

1. Define the analysis and synthesis using the **Hamming window** at 50% overlap.
2. Verify the **reconstruction** of the signal in the time domain using the **overlap-add principle** (OLA).

---

### MATLAB Code for STFT with OLA Reconstruction

```matlab
% Parameters
x = sin(2*pi*0.01*(0:999))'; % Example signal: a sine wave
M = 256; % Window size
R = M/2; % 50% Overlap
w = hamming(M); % Analysis window

% STFT Analysis
num_frames = floor((length(x) - M) / R) + 1;
X = zeros(M, num_frames);
for k = 0:num_frames-1
    x_seg = x(k*R + (1:M)); % Extract frame
    x_seg_windowed = x_seg .* w; % Apply window
    X(:,k+1) = fft(x_seg_windowed); % FFT
end

% ISTFT and Overlap-Add Reconstruction
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
```

**Output**:

- Verify the reconstructed signal visually (should match the original signal).
- Small differences can be analyzed with an error metric like `norm(x - x_reconstructed)`.

---

## **C- Implement Weighted Overlap-Add (WOLA)**

### WOLA: Weighted Overlap-Add introduces **different windows** for analysis and synthesis to improve reconstruction fidelity.

- The analysis window wa[n]w_a[n] is used before the Fourier Transform.
- The synthesis window ws[n]w_s[n] is applied after the Inverse Fourier Transform.
- A commonly-used approach is to use a **squared version** of the window or its complement.

---

### MATLAB Code for WOLA with Unique Windows

```matlab
% Parameters
x = sin(2*pi*0.01*(0:999))'; % Test signal
M = 256; % Window size
R = M/2; % 50% Overlap
w_analysis = hamming(M); % Analysis window
w_synthesis = w_analysis.^2; % Example synthesis window: squared version

% STFT Analysis
num_frames = floor((length(x) - M) / R) + 1;
X = zeros(M, num_frames);
for k = 0:num_frames-1
    x_seg = x(k*R + (1:M)); % Extract frame
    x_seg_windowed = x_seg .* w_analysis; % Apply analysis window
    X(:,k+1) = fft(x_seg_windowed); % FFT
end

% ISTFT and WOLA Reconstruction
x_reconstructed = zeros(length(x), 1); % Initialize
for k = 0:num_frames-1
    x_ifft = ifft(X(:,k+1)); % Inverse FFT
    x_overlap = real(x_ifft) .* w_synthesis; % Apply synthesis window
    x_reconstructed(k*R + (1:M)) = x_reconstructed(k*R + (1:M)) + x_overlap;
end

% Plot results
figure;
plot(x, 'b'); hold on;
plot(x_reconstructed, 'r--');
title('WOLA Reconstruction vs Original');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Samples');
ylabel('Amplitude');
```

### Key Notes:

- The synthesis window can also be optimized using COLA-compatible modifications.
- Compare the reconstruction quality between WOLA and standard OLA.

---

### **Deliverables**:

1. Verification of COLA property for all mentioned windows.
2. Demonstration of OLA-based STFT and ISTFT with a Hamming window.
3. Implementation of WOLA using unique analysis and synthesis windows.

___
To reproduce the figures in your image, we’ll use MATLAB for the following:

1. **Divide a signal into overlapping segments** (using a window with a specific overlap percentage).
2. **Process each segment (e.g., FFT and IFFT)**.
3. **Reconstruct the original signal using the overlap-add method**.

The goal is to simulate the **COLA (Constant-Overlap-Add)** constraint and demonstrate the proper reconstruction.

---

## **COLA (Constant-Overlap-Add)** constraint simulation and reconstruction demonstration

Here's a detailed MATLAB script:

```matlab
% Parameters
fs = 1000;               % Sampling frequency
t = 0:1/fs:5-1/fs;       % Time vector for 5 seconds
x = chirp(t,0,5,250);    % Test signal: Chirp signal from 0 Hz to 250 Hz

M = 512;                 % Window size
overlap = 50;            % Overlap percentage (50% in this case)
R = M * (1 - overlap/100); % Hop size

w = hamming(M);          % Window: Hamming
num_frames = floor((length(x) - M) / R) + 1;

% Plotting setup
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
```

---
## **Steps Explained**:

1. **Generate Input Signal**: A chirp signal is used for visualization.
2. **Windowed Buffers**:
    - Use a **Hamming window** with 50% overlap.
    - Visualize the overlapping windows.
3. **STFT and ISTFT Processing**:
    - Apply FFT to each segment, then an optional signal manipulation step.
    - Perform IFFT to return to the time domain.
4. **Overlap-Add Reconstruction**:
    - Sum all windowed segments back into the time domain to demonstrate the COLA principle.
5. **Plot the Results**:
    - Original signal, windows, processed windows, and the final reconstructed signal.

---
## **Output**:

- **Plot 1**: The original signal, showing the windowed buffers.
- **Plot 2**: Hamming windows with 50% overlap applied to the signal.
- **Plot 3**: The individual processed (FFT/IFFT) windowed segments.
- **Plot 4**: Reconstructed signal using overlap-add.
- **Plot 5**: Comparison of the original and reconstructed signals.

---
## **Verification**:

- The reconstructed signal in **Plot 5** should match the original signal. This validates the COLA property and demonstrates the proper use of the window function with overlap.