function [y_resamples] = polyphase_decimation(x,h,Down_Ratio)
% polyphase decimation filter
% 
% This function performs polyphase decimation, which involves filtering the
% input signal `x` with a given resampling filter `h` and downsampling it by
% the downsampling ratio `Down_Ratio`.
%
% Inputs:
%   x           : Input signal (can be a row or column vector).
%   h           : Resampling filter (FIR filter designed for resampling).
%   Down_Ratio  : The decimation factor, i.e., the factor by which the
%                 signal is downsampled.
%
% Output:
%   y_resamples : The decimated and filtered output signal after polyphase
% filtering.

% First, check if any elements in the filter coefficients `h` are NaN.
% If there are NaNs, truncate the filter to only include the valid
% coefficients up to the first NaN.
for index = 1:length(h)
    if isnan(h(index)) == 1
        h = h(1:index-1); % Truncate the filter at the first NaN
        break;
    end
end

% Reshape the input signal to a row vector (if it's a column vector)
x = x(:)';  % Transpose the input signal to ensure it's a row vector.

% Determine the length of the input signal after downsampling by
% `Down_Ratio`.
% This ensures that the input signal length is a multiple of `Down_Ratio`
% for even downsampling.
l_in = Down_Ratio * floor(length(x) / Down_Ratio);

% Create a downsampled version of `x` where every `Down_Ratio`-th sample is
% kept.
% This will be used to create the polyphase components for filtering.
x_dec(1,:) = [x(1:Down_Ratio:l_in), 0];  % Downsample the signal and append
% a 0 for length consistency.

% For each phase in the polyphase decomposition, downsample the signal
% accordingly.
for index = 2:Down_Ratio
    % For each phase, take every `Down_Ratio`-th sample starting from a
    % different offset.
    x_dec(index,:) = [0, x(Down_Ratio - index + 2:Down_Ratio:l_in)]; 
end

% Polyphase filter decomposition: Split the filter `h` into `Down_Ratio`
% phases.
% This divides the filter into sections that correspond to the different
% downsampling phases.
for index = 1:Down_Ratio
    % Reshape the filter into phases, each phase corresponds to a part of the filter coefficients.
    % The filter is decomposed into blocks of `Down_Ratio` coefficients.
    E(index,:) = reshape(h(index:Down_Ratio:end), 1, length(h(index:Down_Ratio:end))); 
    % `E(index,:)` is the `index`-th phase of the filter, reshaped into a row vector.
end

% Perform the polyphase filtering: For each phase, filter the corresponding
% downsampled signal.
for index_down = 1:Down_Ratio
    % Apply the filter for each phase `index_down` to the corresponding
    % downsampled signal `x_dec(index_down,:)`.
    % This uses the `filter` function to filter the downsampled signal
    % `x_dec(index_down,:)`
    % with the corresponding phase filter `E(index_down,:)`.
    u(index_down,:) = filter(squeeze(E(index_down,:)), 1, x_dec(index_down,:)); 
    % `squeeze` is used to remove any singleton dimensions from the phase filter.
end

% Sum the results of all the polyphase filtering stages to get the final
% decimated and filtered output.
y_resamples = sum(u);  % The output is the sum of all the phases' filtered signals.
