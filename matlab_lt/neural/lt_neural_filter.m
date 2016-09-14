function [datfilt,neuralFiltLow,neuralFiltHi] =lt_neural_filter(dat, frequency_parameters, use_waveclusver)
%% lt 8/17/16 - bandpass filters intan neural data.

% dat, vector
% frequency_parameters, from intan file (or just the fs num,ber, e.g.
% 25000)
% waveclusver=1, then filters same way wave_clus does. else does
% butterworth


% default, use wave_clus version
if ~exist('waveclusver', 'var')
    use_waveclusver=1;
end


neuralFiltLow=300;
neuralFiltHi=3000;

try
    fs=frequency_parameters.amplifier_sample_rate;
catch err
    fs=frequency_parameters;
end

%% ==== make neural filter
if use_waveclusver==0
    %% then do butterworth bandpath
    N=4;
    
    % if neuralFiltLow<min([frequency_parameters.actual_lower_bandwidth, ...
    %         frequency_parameters.actual_dsp_cutoff_frequency])
    %     % make low freq at least as high at acquisition frequency.
    %     neuralFiltLow=min([frequency_parameters.actual_lower_bandwidth, ...
    %         frequency_parameters.actual_dsp_cutoff_frequency]);
    % end
    % neuralFiltHi=frequency_parameters.actual_upper_bandwidth;
    [filt_b,filt_a]=butter(N,[neuralFiltLow*2/fs, neuralFiltHi*2/fs]);
    
    
    % === filter
    datfilt=filtfilt(filt_b,filt_a,dat);
    
elseif use_waveclusver==1
    %% wave_clus version [use this, gets higher SNR for spikes, I think
    
    x=dat;
    
    if exist('ellip','file')                         %Checks for the signal processing toolbox
        [b_detect,a_detect] = ellip(4,0.1,40,[neuralFiltLow 3000]*2/fs);
        [b,a] = ellip(2,0.1,40,[neuralFiltLow 3000]*2/fs);
        %     [z_det,p_det,k_det] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
        %     [z,p,k] = ellip(par.sort_order,0.1,40,[fmin_sort fmax_sort]*2/sr);
        %
        %     [SOS,G] = zp2sos(z,p,k);
        %     [SOS_det,G_det] = zp2sos(z_det,p_det,k_det);
        if exist('FiltFiltM','file')
            xf_detect = FiltFiltM(b_detect, a_detect, x);
            xf = FiltFiltM(b, a, x);
        else
            xf_detect = filtfilt(b_detect, a_detect, x);
            xf = filtfilt(b, a, x);
            %         xf_detect = filtfilt(SOS_det, G_det, x);
            %         xf = filtfilt(SOS,G, x);
            
        end
        
    else
        xf = fix_filter(x);                   %Does a bandpass filtering between [300 3000] without the toolbox.
        xf_detect = xf;
    end
    
    datfilt=xf;
end

