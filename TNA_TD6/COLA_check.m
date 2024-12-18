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