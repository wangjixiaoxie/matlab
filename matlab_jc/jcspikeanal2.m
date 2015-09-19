% jcspikeanal('b',0,0,'combined_data...','o',0.0400,premotor_window);
% 'b' is the note or context to search for
function [MAT_spiketimes,MAT_nspikes,pitch]=jcspikeanal2(search_for,skip_first,skip_last,loadfilename,symbol,t_assay,premotor_window)
eval(sprintf('load %s',loadfilename))

premotor_estimate=premotor_window;

spiketimes=round(spiketimes*1000);   % in milliseconds
spiketimes=unique(spiketimes);  % if more than one spike in a bin, call it one spike

sp=zeros(1,length(spiketimes));

spiketimes=nonzeros(spiketimes);

sp(spiketimes)=1000;

onsets=unique(round(onsets));
offsets=unique(round(offsets));


% load for pitch information




% find notes
for x=1:(length(labels)-length(search_for)+1)
    compare=labels(x:(x+length(search_for)-1));
    true(x)=strcmp(search_for,compare);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get spike number and syl duration
id_true=find(true);                                             % id of starts of sequence
%id_true,search_for,return
disp(['Found ' num2str(length(id_true))  ' sequences'])
MAT_duration=zeros(length(id_true),length(search_for));
MAT_int_duration=MAT_duration;
MAT_nspikes=MAT_duration;
MAT_nspikes_int=MAT_duration;
for x=1:length(id_true)
    for y=1:length(search_for)

        % for calculating syllable durations
        syllable_onsets=onsets(id_true(x)+y-1);
        syllable_offsets=offsets(id_true(x)+y-1);
        if x<length(id_true)
            syllable_onsets_nextsyl=onsets(id_true(x)+y);
        else
            syllable_onsets_nextsyl=offsets(id_true(x)+y-1);
        end

        % this is where window for spike counting is set
        %disp(['At window setting - [pre_onset post_onset] = ' num2str([pre_onset post_onset])])
        % OLD  - this was wrong - pre_ and post_onset should be relative to
        % time of quantification, not to onset of syllable
        %         on=onsets(id_true(x)+y-1)-pre_onset;
        %         off=onsets(id_true(x)+y-1)+post_onset;

        if t_assay<1
            pre_onset=-1*(t_assay*1000-premotor_estimate(2));
            post_onset=t_assay*1000-premotor_estimate(1);
        else   % if t_assay>1, it is the percent time through the syllable.
            t_assay_with_pct=.01*t_assay*[offsets(id_true(x)+y-1)-onsets(id_true(x)+y-1)];
            pre_onset=-1*(t_assay_with_pct-premotor_estimate(2));
            post_onset=t_assay_with_pct-premotor_estimate(1);
            % %assuming that avg duration is 50ms.
            %             pre_onset=-1*(.01*t_assay*50-premotor_estimate(2));
            %             post_onset=.01*t_assay*50-premotor_estimate(1);
        end

        % this is where window for spike counting is set
        on=onsets(id_true(x)+y-1)-pre_onset;
        off=onsets(id_true(x)+y-1)+post_onset;

        if numel(onsets)>=(id_true(x)+y-1)
            MAT_nspikes(x,y)=length(find(spiketimes>on & spiketimes<off));   % number of spikes between on and off;
            MAT_spiketimes{x,y}=spiketimes(find(spiketimes>on & spiketimes<off))-onsets(id_true(x)+y-1)-t_assay*1000;
            % old measures when on and off corresponded to onset and offset
            % of syl
            MAT_duration(x,y)=(syllable_offsets-syllable_onsets);
            MAT_int_duration(x,y)=(syllable_onsets_nextsyl-syllable_offsets);

        end
    end
    MAT_amplitude(x,:)=peak_amplitudes(id_true(x):(id_true(x)+length(search_for)-1));
    MAT_amplitude_summed(x,:)=summed_amplitudes(id_true(x):(id_true(x)+length(search_for)-1));
    if exist('amplitude_at_pitchquant')
        MAT_amplitude_at_pitchquant(x,:)=amplitude_at_pitchquant(id_true(x):(id_true(x)+length(search_for)-1));
    end
    MAT_amplitude_16msec_around_pitchquant(x,:)=amplitude_16msec_around_pitchquant(id_true(x):(id_true(x)+length(search_for)-1));
    MAT_logamp(x,:)=log10(MAT_amplitude_16msec_around_pitchquant(x,:));
end

