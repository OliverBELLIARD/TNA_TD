function [y_resamples, y_r] = polyphase_fractional(x,h,Up_Ratio,Down_Ratio,output_delay)
% Fractional resampling using fractional polyphase filter
%
%   [y_resamples, y_r] = polyphase_fractional(x,h,Up_Ratio,Down_Ratio,output_delay)
%
% Inputs:
%   x             : Input signal (could be a row or column vector).
%   h             : Resampling filter coefficients (FIR filter).
%   Up_Ratio      : Upsampling ratio (factor by which the signal is upsampled).
%   Down_Ratio    : Downsampling ratio (factor by which the signal is downsampled).
%   output_delay  : Vector of delays for each phase in the filter.
%
% Outputs:
%   y_resamples   : The resampled output signal after fractional resampling
%                   (both upsampling and downsampling).
%   y_r           : The upsampled signal after filtering.

% Truncate the filter `h` if it contains NaN values.
for index = 1:length(h)
    if isnan(h(index)) == 1
        h = h(1:index-1); % If NaN is found, truncate the filter up to the first NaN element.
        break;
    end
end

% Reshape the input signal `x` to be a row vector (if it's not already).
x = x(:)';  % Convert `x` to a row vector.

% Compute the length of the input signal after adjusting for the
% downsampling ratio.
l_in = Down_Ratio * floor(length(x)/Down_Ratio);  % Ensures that the signal length is a multiple of `Down_Ratio`.

% Initialize the downsampled signal matrix `x_dec` to store the downsampled
% signal in multiple phases.
x_dec(1,:) = [x(1:Down_Ratio:l_in), 0];  % Downsample the signal by `Down_Ratio` and append 0 for length consistency.

% For each of the remaining phases, downsample the signal with different offsets.
for index = 2:Down_Ratio
    % For each phase, start the downsampling from a different offset.
    x_dec(index,:) = [0, x(Down_Ratio - index + 2:Down_Ratio:l_in)];
end

% Resampler decomposition: Decompose the filter `h` into polyphase components.
for index = 1:Down_Ratio
    % Reshape the filter into `Up_Ratio` phases. Each phase corresponds to
    % a section of the filter.
    % The filter coefficients are spread across `Up_Ratio` rows for each phase.
    E(index,:,:) = reshape(h(index:Down_Ratio:end), Up_Ratio, length(h(index:Down_Ratio:end)) / Up_Ratio); 
end

% Perform polyphase filtering and upsampling:
% For each phase of the downsampled signal, apply the corresponding filter
% phase and upsample.
for index_down = 1:Down_Ratio
    for index_up = 1:Up_Ratio
        % Filter the downsampled signal and upsample it by multiplying by `Up_Ratio`.
        u(index_down,index_up,:) = Up_Ratio * filter(squeeze(E(index_down,index_up,:)), 1, x_dec(index_down,:));
    end
end

% Determine the size of the filtered output `u` for further processing.
[l, vec_len] = size(u);  % `l` is the number of downsampling phases, `vec_len` is the length of each filtered signal.

% Compute the maximum delay across all phases.
max_output_delay = max(output_delay);  % Find the maximum delay in the output delays.

% Initialize a matrix `yk` to store the delayed filtered output for each phase.
yk = zeros(Down_Ratio, (max_output_delay + vec_len));  % Allocate memory for delayed results.

% Apply the delays for each phase and store the filtered signals in `yk`.
for index_down = 1:Down_Ratio
    bound_low = output_delay(index_down) + 1;  % Lower bound for the delayed output.
    bound_high = bound_low + vec_len - 1;      % Upper bound for the delayed output.
    
    % Extract the filtered output for the current downsampling phase and flatten it.
    y_temp = squeeze(u(index_down,:,:));
    y_col = y_temp(:)';  % Flatten the phase output into a row vector.
    
    % Insert the delayed filtered output into the corresponding positions in `yk`.
    yk(index_down, bound_low:bound_high) = y_col;
end

% Sum the filtered outputs across all downsampling phases to get the final
% resampled signal.
if Down_Ratio ~= 1 
    y_resamples = sum(yk);  % If Down_Ratio is not 1, sum the phases to get the final output.
else
    y_resamples = yk;  % If Down_Ratio is 1, no downsampling is applied, return `yk` as is.
end

% Perform upsampling of the input signal `x`.
n = 1:length(x) * Up_Ratio;  % Create an index for the upsampled signal.
x_up = zeros(1, length(x) * Up_Ratio);  % Initialize the upsampled signal.
x_up(1:Up_Ratio:end) = x;  % Insert the input signal `x` into every `Up_Ratio`-th position.

% Apply the resampling filter `h` to the upsampled signal.
y_up = Up_Ratio * filter(h, 1, x_up);  % Filter the upsampled signal and scale by `Up_Ratio`.

% Downsample the filtered signal by `Down_Ratio` to produce the final
% fractional resampled signal.
y_r = y_up(1:Down_Ratio:end);  % Downsample by taking every `Down_Ratio`-th sample.
end
