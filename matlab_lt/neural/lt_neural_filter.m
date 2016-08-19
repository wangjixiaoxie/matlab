function [datfilt,neuralFiltLow,neuralFiltHi] =lt_neural_filter(dat, frequency_parameters)
%% lt 8/17/16 - bandpass filters intan neural data.

% dat, vector
% frequency_parameters, from intan file


neuralFiltLow=500;
fs=frequency_parameters.amplifier_sample_rate;

%% ==== make neural filter
% butterworth bandpath
N=4;

if neuralFiltLow<min([frequency_parameters.actual_lower_bandwidth, ...
        frequency_parameters.actual_dsp_cutoff_frequency])
    % make low freq at least as high at acquisition frequency.
    neuralFiltLow=min([frequency_parameters.actual_lower_bandwidth, ...
        frequency_parameters.actual_dsp_cutoff_frequency]);
end

neuralFiltHi=frequency_parameters.actual_upper_bandwidth;

[filt_b,filt_a]=butter(N,[neuralFiltLow*2/fs, neuralFiltHi*2/fs]);


% === filter

            datfilt=filtfilt(filt_b,filt_a,dat);
