if 0
    canonical_onsets=onsets_all_syls{1};
    canonical_offsets=offsets_all_syls{1};
elseif 1
    canonical_onsets=mean(onsets_all_syls_MAT,1);
    canonical_offsets=mean(offsets_all_syls_MAT,1);
    if 0
        disp('Using mean onset and offset time')
        eval(sprintf('save canonical_vec_3_17_2006_1409-1544_%s canonical_onsets canonical_offsets',search_for))
    end
else
    load_twarp_filename=['canonical_vec_3_17_2006_1409-1544_' search_for];
    disp(['Loading time warp file ' load_twarp_filename])
    eval(sprintf('load canonical_vec_3_17_2006_1409-1544_%s ',search_for))
end

ct=1;
for x=1:length(canonical_onsets);
    canonical_vec(ct)=canonical_onsets(x);ct=ct+1;  % canonical vec is on(1) off(1) on(2) off(2) etc
    canonical_vec(ct)=canonical_offsets(x);ct=ct+1;
end

spiketimes_motif=spiketimes(find(spiketimes>canonical_vec(1) & spiketimes<canonical_vec(end))); %spiketimes between these

canonical_vec=canonical_vec-canonical_vec(1);   % set first time to zero

n=1:length(onsets_all_syls);
%clear global stretch_pct
%global stretch_pct

for yy=1:length(n)%n_trials_to_use%1:length(onsets_all_syls) % non-canonical - warp these
    y=n(yy);
    start_motif=onsets_all_syls{y}(1);   %first onset
    end_motif=offsets_all_syls{y}(end);   % last offset
    
    spiketimes_motif=spiketimes(find(spiketimes>start_motif & spiketimes<end_motif)); %spiketimes between these
    spiketimes_motif=spiketimes_motif-start_motif;

    ct=1;
    for x=1:length(canonical_onsets);
        vec_to_be_warped(ct)=onsets_all_syls{y}(x);ct=ct+1;  % canonical vec is on(1) off(1) on(2) off(2) etc
        vec_to_be_warped(ct)=offsets_all_syls{y}(x);ct=ct+1;
    end
    vec_to_be_warped=vec_to_be_warped-vec_to_be_warped(1);
    stretch_pct(yy,:)=1-diff(vec_to_be_warped)./diff(canonical_vec);
    spiketimes_warped{yy}= interp1(vec_to_be_warped,canonical_vec,spiketimes_motif,'linear');
% the unwarped spiketimes of the FIRST instance of the motif
% this is used to put dots on spikes in the example neural trace
    if yy==1
        first_spiketrain_unwarped=spiketimes_motif;
    end
    if(sum(isnan(spiketimes_warped{yy})))
        disp('WARNING - NAN')
    end
end
