%% tested (300-3000hz) performance of butter (old code) and elipse from wave_clus
%%  for visualization I am still using butter, but for wave_clus I think other code is better
% better peaks with similar amount of noise.
%%  do not run this code. just for my notes.
%%  also tried using high top freq (like 14000) but not much difference. lower (like 3000) is 
% probably better

%%  butter version.

% butterworth bandpath
N=4;

if neuralFiltLow<min([frequency_parameters.actual_lower_bandwidth, ...
        frequency_parameters.actual_dsp_cutoff_frequency])
    % make low freq at least as high at acquisition frequency.
    neuralFiltLow=min([frequency_parameters.actual_lower_bandwidth, ...
        frequency_parameters.actual_dsp_cutoff_frequency]);
end

neuralFiltHi=frequency_parameters.actual_upper_bandwidth;
neuralFiltHi=3000;

[filt_b,filt_a]=butter(N,[neuralFiltLow*2/fs, neuralFiltHi*2/fs]);

            dat=filtfilt(filt_b,filt_a,amplifier_data(i, :));

            
           %% ==== OR DO WAVE_CLUS VERSION
           % version 1 and version 2 are the same. version 1 is not using
           % signal toolbox. not sure if is [300 3000] or [400 3000];
            x=amplifier_data(i, :);
            a = [1.0000 -2.3930  2.0859 -0.9413 0.2502];
            b = [0.1966 -0.0167 -0.3598 -0.0167 0.1966];
            
            x = x(:);
            len = size(x,1);
            b = b(:).';
            a = a(:).';
            nb = length(b);
            na = length(a);
            nfilt = max(nb,na);
            
            nfact = 3*(nfilt-1);  % length of edge transients
            
            rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
            cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
            data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
            sp = sparse(rows,cols,data);
            zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
            
            y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];
            
            if exist('FilterM','file')
                y = FilterM(b,a,y,[zi*y(1)]);
            else
                y = filter(b,a,y,[zi*y(1)]);
            end
            y = y(length(y):-1:1);
            
            %second filter, in the other way
            if exist('FilterM','file')
                y = FilterM(b,a,y,[zi*y(1)]);
            else
                y = filter(b,a,y,[zi*y(1)]);
            end
            y = y(length(y):-1:1);
            
            y([1:nfact len+nfact+(1:nfact)]) = [];
            
            y = y.';
            
            
            % =============================
            % another wave_clus version [same as above]
            if exist('ellip','file')                         %Checks for the signal processing toolbox
                [b_detect,a_detect] = ellip(4,0.1,40,[300 3000]*2/fs);
                [b,a] = ellip(2,0.1,40,[300 3000]*2/fs);
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
            