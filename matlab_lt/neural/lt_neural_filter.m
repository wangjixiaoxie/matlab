function [datfilt,neuralFiltLow,neuralFiltHi] =lt_neural_filter(dat, frequency_parameters)
%% lt 8/17/16 - bandpass filters intan neural data.

% dat, vector
% frequency_parameters, from intan file (or just the fs num,ber, e.g.
% 25000)


neuralFiltLow=300;

try
fs=frequency_parameters.amplifier_sample_rate;
catch err
    fs=frequency_parameters;
end

%% ==== make neural filter
% butterworth bandpath
N=4;

% if neuralFiltLow<min([frequency_parameters.actual_lower_bandwidth, ...
%         frequency_parameters.actual_dsp_cutoff_frequency])
%     % make low freq at least as high at acquisition frequency.
%     neuralFiltLow=min([frequency_parameters.actual_lower_bandwidth, ...
%         frequency_parameters.actual_dsp_cutoff_frequency]);
% end
neuralFiltLow=300
% neuralFiltHi=frequency_parameters.actual_upper_bandwidth;
neuralFiltHi=3000;
[filt_b,filt_a]=butter(N,[neuralFiltLow*2/fs, neuralFiltHi*2/fs]);


% === filter
datfilt=filtfilt(filt_b,filt_a,dat);
