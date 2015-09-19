% wav_file_name is xxx.wav, the name of the wav file containing bos.  will
% load wav_file_name.not.mat.  the .not.mat file was generated with
% BOS_make_combined_data.m
function foo(search_for,loadfile,wav_file_name)

quantify_amp_and_pitch=1;

spike_color='b';
if ~isstr(loadfile)
    ch_with_data=loadfile;
else
    tmp=strfind(loadfile,'CH');
    ch_with_data=str2num(loadfile(tmp+2));
end
sort_by_pitch=0;
mark_intervals_and_pitchquant=1;    % mark pitchquant doesnt work yet
pitch_of_syl_number=1;
use_preceding_interval=1;
plot_spect_and_neural_ex=1;
lines_connecting_songs=0;

neural_available=1;
plot_spect_and_neural_ex=plot_spect_and_neural_ex*neural_available;
if strfind(pwd,'pu44w52')
    birdname='pu44w52';
    song_channel{1}='obs0';
    disp(' ');        disp(' This is for bird pu44w52 - using song channel obs0 ');        disp(' '); pause(1)

elseif strfind(pwd,'pk43gr64')
    birdname='pk43gr64';
    song_channel{1}='obs2';
    disp(' ');        disp(' This is for bird pk43gr64 - using song channel obs2 ');        disp(' '); pause(1)
elseif strfind(pwd,'g26g23')
    birdname='g26g23';
    song_channel{1}='obs0r';
    disp(' ');        disp(' This is for bird g26g23 - using song channel obs0r ');        disp(' '); pause(1)
elseif strfind(pwd,'pu24w39')
    birdname='pu24w39';
    song_channel{1}='obs0';
    disp(' ');        disp(' This is for bird pu24w39');        disp(' '); pause(1)
elseif strfind(pwd,'r75g59')
    birdname='r75g59';
    disp(' ');        disp(' This is for bird r75g59');        disp(' '); pause(1)
    song_channel{1}='obs0r';
elseif strfind(pwd,'pu5b3')
    birdname='pu5b3';
    disp(' ');        disp(' This is for bird pu5b3');        disp(' '); pause(1)
    song_channel{1}='obs0r';
elseif strfind(pwd,'w48o28')
    birdname='w48o28';
    disp(' ');        disp(' This is for bird w48o28');        disp(' '); pause(1)
    song_channel{1}='obs0r';
elseif strfind(pwd,'pu54b57')
    birdname='pu54b57';
    disp(' ');        disp(' This is for bird pu54b57');        disp(' '); pause(1)
    song_channel{1}='obs0r';
elseif strfind(pwd,'pu53b58')
    birdname='pu53b58';
    disp(' ');        disp(' This is for bird pu53b58');        disp(' '); pause(1)
    song_channel{1}='obs0r';
elseif strfind(pwd,'pu28b75')
    birdname='pu28b75';
    disp(' ');        disp(' This is for bird pu28b75');        disp(' '); pause(1)
    song_channel{1}='obs0r';
elseif strfind(pwd,'g18g8')
    birdname='g18g8';
    disp(' ');        disp(' This is for bird g18g8');        disp(' '); pause(1)
    song_channel{1}='obs0';
elseif strfind(pwd,'ZFpi90w37')
    birdname='ZFpi90w37';
    disp(' ');        disp(' This is for bird ZFpi90w37');        disp(' '); pause(1)
    song_channel{1}='obs0';
elseif strfind(pwd,'g44o23')
    birdname='g44o23';
    disp(' ');        disp(' This is for bird g44o23');        disp(' '); pause(1)
    song_channel{1}='obs0';
elseif strfind(pwd,'pu26y2')
    birdname='pu26y2';
    disp(' ');        disp(' This is for bird pu26y2');        disp(' '); pause(1)
    song_channel{1}='obs0';
elseif strfind(pwd,'r12r11')
    birdname='r12r11';
    disp(' ');        disp(' This is for bird r12r11');        disp(' '); pause(1)
    song_channel{1}='obs0';
elseif strfind(pwd,'B53O71')    % Mel's bird
    birdname='B53O71';
    disp(' ');        disp(' This is for bird B53O71');        disp(' '); pause(1)
    song_channel{1}='obs4';
elseif strfind(pwd,'O55Pu53')    % Mel's bird
    birdname='O55Pu53';
    disp(' ');        disp(' This is for bird O55Pu53');        disp(' '); pause(1)
    song_channel{1}='obs4';
elseif strfind(pwd,'G26G23')    % Mel's bird
    birdname='G26G23';
    disp(' ');        disp(' This is for bird G26G23');        disp(' '); pause(1)
    song_channel{1}='obs4';
elseif strfind(pwd,'G45G46')    % Mel's bird
    birdname='G45G46';
    disp(' ');        disp(' This is for bird G45G46');        disp(' '); pause(1)
    song_channel{1}='obs2';
elseif strfind(pwd,'Pu55Pu22')    % Mel's bird
    birdname='Pu55Pu22';
    disp(' ');        disp(' This is for bird Pu55Pu22');        disp(' '); pause(1)
    song_channel{1}='obs4';
end
if  isnumeric(loadfile)




    if isunix
        !rm batchfile
        !ls -1 *bin > batchfile
    else
        !dir /B **.cbin > batchfile
    end
    fid=fopen('batchfile','r');

    song=[];
    labels_combo=[];
    onsets_combo=[];
    offsets_combo=[];
    spiketimes_combo=[];
    old_time=0;
    peak_amp_combo=[];
    summed_amp_combo=[];
    amplitude_at_pitchquant_combo=[];
    amplitude_16msec_around_pitchquant_combo=[];
    peak_labelvec_combo=[];
    peak_pinterp_labelvec_combo=[];
    pct_error=[];
pct_ISIs_leq_1ms=[];
   recommended_TH_combo=[];

    first_ch=1;
    ct=1;
    while 1
        fn=fgetl(fid);
        if (~ischar(fn));break;end
        disp(['Loading ' fn '.not.mat and .neuralnot_CH' num2str(ch_with_data) '.mat'])
        clear Fs onsets fname peak_labelvec labels peak_pinterp_labelvec  min_dur sm_win min_int  threshold  offsets
        clear peak_labelvec amplitude_at_pitchquant 
        clear peak_amp summed_amp  amplitude_16msec_around_pitchquant
        eval(sprintf('load %s.not.mat',fn))
        
        if neural_available
            eval(sprintf('load %s.neuralnot_CH%s.mat',fn,num2str(ch_with_data)))
            %            eval(sprintf('load %s.neuralnot.mat',fn))
        else
            Data.total_n_samples=0;
            Data.channels=1;
            Data.SNR=0;
            Data.Fs=32000;
            Data.spiketimes{1}=1;
        end
        if quantify_amp_and_pitch
        if length(peak_pinterp_labelvec)~=length(labels);error('lengths dont agree');end
            if strcmp(birdname,'pu44w52')
                id_correct_syl_b=find(labels=='b' & peak_pinterp_labelvec<5000);
                if ~isempty(id_correct_syl_b)
                    labels(id_correct_syl_b)='x';
                    disp(' ')
                    for xyz=1:length(id_correct_syl_b)
                        disp('Correcting syllable "b" by renaming as syllable "x" if pitch<5000')
                    end
                    disp(' ')
                end
            elseif strcmp(birdname,'g26g23')    % set syl 'b' pitch = 0 - for MEL's gg bird
                id_b=find(labels=='b');
                peak_pinterp_labelvec(id_b)=0;
                disp(' ');disp('Setting pitch of syllable "b" = 0 - DO ONLY FOR MELS GG BIRD');disp(' ')
            else
                disp(' ');disp('No syllable corrections being implemented');disp(' ')
            end
        end

        if first_ch & isempty(input(['Use recommended threshold to get spiketimes? (return for yes, "n" for no)'],'s'));
            if ~exist('use_recommended_TH');use_recommended_TH=1;end
            first_ch=0;
        else
            if ~exist('use_recommended_TH');use_recommended_TH=0;end
            if Data.n_units>1 & first_ch
                if ~exist('loadfile') | ~isstr(loadfile)
                    use_unit=input(['Spike data available for ' num2str(Data.n_units) ' units.  use which? (larger number = larger unit), or return to use TH_recommended ']);
                    use_unit_entry_str=num2str(use_unit);
                end
                first_ch=0;
            elseif first_ch
                %            use_unit=Data.channels;
                use_unit=1;
                first_ch=0;
            end
        end
        
if quantify_amp_and_pitch
        for x=1:length(song_channel)
            disp(['Loading song to compute peak amplitude'])
            [song, Fs_sound]=soundin('.', fn,song_channel{1});
        end

        % compute peak amplitude
        for x=1:length(onsets)
            onsyl=round(onsets(x)*32);  % in sample number (onsets is in msec)
            if onsyl==0;onsyl=1;disp('Setting onsyl=1 in neural_by_syl_sequence');end
            offsyl=round(offsets(x)*32);
            if offsyl>length(song);offsyl=length(song);end
            sound_clip=song(onsyl:offsyl);
            [smooth,spec,t,f]=evsmooth(sound_clip,32000);
            %            figure;plot(smooth)
            peak_amp(x)=max(smooth);
            summed_amp(x)=sum(smooth);
            amplitude_16msec_around_pitchquant(x)=compute_summed_amp_window(smooth,birdname,labels(x));
        end
        peak_amp_combo=[peak_amp_combo peak_amp];
        summed_amp_combo=[summed_amp_combo summed_amp];
        amplitude_16msec_around_pitchquant_combo=[ amplitude_16msec_around_pitchquant_combo  amplitude_16msec_around_pitchquant];
        %    [neural_tmp, Fs]=soundin('.', fn,['obs' num2str(use_unit)]);
        if Fs==3  | Fs_sound==3;Fs=32000;Fs_sound=32000;end   % I still dont understant this bug
else
    disp('NOT loading pitches or computing peak amplitude')
end

        labels_combo=[labels_combo labels];
        L_labels(ct)=length(labels);
        fname_arr{ct}=fn;
        onsets_combo=[onsets_combo; onsets+sum(old_time)*1000];
        offsets_combo=[offsets_combo; offsets+sum(old_time)*1000];

        if quantify_amp_and_pitch
        peak_labelvec_combo=[peak_labelvec_combo peak_labelvec];
        peak_pinterp_labelvec_combo=[peak_pinterp_labelvec_combo peak_pinterp_labelvec];
        amplitude_at_pitchquant_combo=[amplitude_at_pitchquant_combo amplitude_at_pitchquant];
        end
        if use_recommended_TH
            disp('Using recommended threshold for spike times')
            spiketimes_combo=[spiketimes_combo Data.spiketimes_from_recommended_TH + sum(old_time)];
        else
            disp('Using cluster for spike times')
            spiketimes_combo=[spiketimes_combo Data.spiketimes{use_unit+1}+sum(old_time)];
        end

        recommended_TH_combo=[recommended_TH_combo Data.recommended_TH];

        total_samples_cat(ct)=Data.total_n_samples;

        if use_recommended_TH
            % Data.pct_error is of length =  number of clusters 
            % Data.pct_ISIs_leq_1ms is of length =  number of clusters + 1,
            % where the last entry is pct_ISIs_leq_1ms for the thresholded
            % spiketimes
            pct_error=[pct_error; Data.pct_error(end)];
            pct_ISIs_leq_1ms=[pct_ISIs_leq_1ms; Data.pct_ISIs_leq_1ms(end)];
        else
            pct_error=[pct_error; Data.pct_error(use_unit+1)];
            pct_ISIs_leq_1ms=[pct_ISIs_leq_1ms; Data.pct_ISIs_leq_1ms(use_unit+1)];
        end

        %        old_time(ct)=Data.total_n_samples/Data.Fs;
        old_time(ct)=Data.total_n_samples/Data.Fs;

        ct=ct+1;
    end
    length_song=sum(total_samples_cat);
    fclose(fid);
    labels=labels_combo;
    onsets=onsets_combo;
    offsets=offsets_combo;
    spiketimes=spiketimes_combo;
recommended_TH=recommended_TH_combo;    
    % below fields will be empty if quantify_amp_and_pitch=0;
    peak_amplitudes=peak_amp_combo;
    summed_amplitudes=summed_amp_combo;
    amplitude_at_pitchquant=amplitude_at_pitchquant_combo;
    amplitude_16msec_around_pitchquant= amplitude_16msec_around_pitchquant_combo;
    peak_labelvec=peak_labelvec_combo;
    peak_pinterp_labelvec=peak_pinterp_labelvec_combo;

    if use_recommended_TH
        extra_str='_TH_recommended';
        use_unit=[];
    else
        if Data.n_units>1
            % largest unit will be #1, smaller with be #2, etc
            unit_number=find(fliplr([1:Data.n_units])==use_unit);
            extra_str=['_UNIT_' num2str(unit_number)];
        else
            extra_str=[];
        end
    end

    eval(sprintf('save combined_data_PCA_CH%s%s use_unit labels onsets offsets spiketimes  L_labels fname_arr Fs  length_song peak_amplitudes summed_amplitudes peak_labelvec peak_pinterp_labelvec pct_error amplitude_at_pitchquant  amplitude_16msec_around_pitchquant pct_ISIs_leq_1ms recommended_TH birdname',num2str(loadfile),extra_str))
    
    if ~neural_available
        return
    end

elseif isstr(loadfile)
    eval(sprintf('load %s',loadfile))
else
    error('here')
end
tminmax_spect=[0 length_song/Fs];


%tmax=length(song)/Fs;
tmax=length_song/Fs;
%time=(1/Fs):(1/Fs):tmax;

search_num=(search_for);
labels_num=(labels); %converts from ascii to numbers



for x=1:(length(labels)-length(search_for)+1)
    compare=labels_num(x:(x+length(search_for)-1));
    true(x)=strcmp(search_for,compare);
end

id_true=find(true);
onsets_true=onsets(id_true)/1000;
offsets_true=offsets(id_true+length(search_for)-1)/1000;
timing_criterion=0;
if timing_criterion    % change id_true to reflect some timing criterion
    for x=1:length(id_true)
        % rows of next_interval_true = length(id_true)
        % columns of id_true = length(search_for)
        % row x is the interval between syllable c and c+1 on iteration x of
        % sequence, when sequence consists of c syllables
        id=id_true(x):(id_true(x)+length(search_for)-1);
        next_interval_true(x,:)=(onsets(id+1)-offsets(id))/1000;
        prev_interval_true(x,:)=(onsets(id)-offsets(id-1))/1000;

        duration(x,:)=(offsets(id)-onsets(id))/1000;
        s_pitch(x,:)=peak_pinterp_labelvec(id);

    end
    % to see dist
    % USE END-1 for last before final i
    %ID_interval=length(search_for)-1;
    ID_interval=1;

    next_int=next_interval_true(:,ID_interval);
    median_interval=median(next_int);

    prev_int=prev_interval_true(:,ID_interval);
    median_interval_prev=median(prev_int);

    duration_syl=duration(:,ID_interval);
    median_duration=median(duration_syl);

    pitch_syl=s_pitch(:,ID_interval);
    median_pitch=median(pitch_syl);

    %next_int_less_than_p1=next_int(find(next_int<.1));
    %next_int=next_int_less_than_p1;disp('only considering intervals less than .1 sec')

    data_to_segretate=next_int;disp('Separating data based on "next_int"')
    threshold=median_interval;disp('Using median of "next_int" as threshold')

    %data_to_segretate=prev_int;disp('Separating data based on "prev_int"')
    %threshold=median_interval_prev;disp('Using median of "prev_int" as threshold')

    %data_to_segretate=duration_syl;disp(' ');disp('Separating data based on "duration"')
    %threshold=median_duration;disp('Using median of "duration" as threshold');disp(' ')

    %data_to_segretate=pitch_syl;disp(' ');disp('Separating data based on "pitch"')
    %threshold=median_pitch;disp('Using median of "pitch" as threshold');disp(' ')

    disp(['Median following duration is ' num2str(threshold)])

    %    threshold=.1;
    greater_than=1;
    use_quartile=0;

    if greater_than
        if use_quartile
            half=data_to_segretate(find(data_to_segretate>threshold));
            median_interval=median(half);
            threshold=median_interval;
        end
        ineqal_str=['Interval greater than ' num2str(threshold)];
    else
        if use_quartile
            half=data_to_segretate(find(data_to_segretate<threshold));
            median_interval=median(half);
            threshold=median_interval;
        end
        ineqal_str=['Interval less than ' num2str(threshold)];
    end
    disp(' ');        disp(ineqal_str);disp(' ');

    %     first=next_interval_true(:,1);
    %     last=next_interval_true(:,length(search_for)-1);
    %     plot(first,last,'.');return

    %       hist(next_interval_true(:,ID_interval),100),
    %       hold on;yl=get(gca,'ylim');plot([1 1]*threshold,yl,'k:','linew',2)
    %       xlabel('Interval length (sec)')
    %       ylabel('# occurrences','fontsize',16,'fontweight','bold')
    %       set(gca,'fontsize',16,'fontweight','bold')
    %       return

    id_true_2=[];
    for x=1:length(id_true)
        %        if next_interval_true(x,ID_interval)>.0657 % for 1544
        if greater_than
            %            if next_interval_true(x,ID_interval)>threshold% & next_interval_true(x,ID_interval)<.1
            if data_to_segretate(x)>threshold% & next_interval_true(x,ID_interval)<.1
                id_true_2=[id_true_2 id_true(x)];
            end
        else
            %            if next_interval_true(x,ID_interval)<threshold %& next_interval_true(x,ID_interval)<.1
            if data_to_segretate(x)<threshold %& next_interval_true(x,ID_interval)<.1
                id_true_2=[id_true_2 id_true(x)];
            end
        end
    end
    disp([' Of ' num2str(length(id_true)) ' correct sequences, keeping ' num2str(length(id_true_2)) ' based on timing criterion'])

    id_true=id_true_2;
else
    ineqal_str=[];
end

if sort_by_pitch
    pitch_syl_id_true=peak_pinterp_labelvec(id_true+pitch_of_syl_number-1);
    original_first_id=id_true(1);
    [tmp,pitch_order]=sort(pitch_syl_id_true);
    id_true=id_true(pitch_order);
    new_id_original_first_id=find(id_true==original_first_id);
else
    new_id_original_first_id=1;
end

%
% if 1% add 60 msec at beginning
%     onsets_true=[onsets(id_true(1))/1000-.060;  onsets(id_true)/1000];
%     offsets_true=[onsets(id_true(1))/1000; offsets(id_true+length(search_for)-1)/1000];
%     id_true=[id_true(1)-1 id_true];
% end
%[onsets_true offsets_true]

labels_sum= cumsum(L_labels);
for x=1:(length(id_true))
    %for x=1:(length(id_true)-1)
    sequence_from_file_number(x)=min(find(id_true(x)<=labels_sum));
end

if use_preceding_interval% use inter-syl interval BEFORE target seq
    for x=1:length(id_true);
        %         onsets_all_syls{x}=onsets(id_true(x):(id_true(x)+length(search_for)-1))/1000;
        %         offsets_all_syls{x}=offsets(id_true(x):(id_true(x)+length(search_for)-1))/1000;
        onsets_all_syls{x}=onsets(id_true(x):(id_true(x)+length(search_for)-1))/1000;
        offsets_all_syls{x}=offsets(id_true(x)-1:(id_true(x)+length(search_for)-1))/1000;
        onsets_all_syls_MAT(x,:)= onsets_all_syls{x}-onsets_all_syls{x}(1);
        offsets_all_syls_MAT(x,:)= offsets_all_syls{x}-onsets_all_syls{x}(1);
    end
else    % old way
    for x=1:length(id_true);
        onsets_all_syls{x}=onsets(id_true(x):(id_true(x)+length(search_for)-1))/1000;
        offsets_all_syls{x}=offsets(id_true(x):(id_true(x)+length(search_for)-1))/1000;
        onsets_all_syls_MAT(x,:)= onsets_all_syls{x}-onsets_all_syls{x}(1);
        offsets_all_syls_MAT(x,:)= offsets_all_syls{x}-onsets_all_syls{x}(1);
    end
end

if     use_preceding_interval
    timewarping_extra_interval
else
    timewarping
end

nsquare=ceil((length(id_true)+1)^.5);
tbin_for_hist=.005;   % sec
%tbin_for_hist=.01   % sec
%    fs_for_hist=10;
MAT=sparse(length(id_true),ceil(max(canonical_vec)/tbin_for_hist));

%    ax(1)=subplot(8,1,1:4);hold on;set(gca,'fontsize',12,'fontname','arial')
ax(1)=subplot(8,1,1:2);hold on;set(gca,'fontsize',12,'fontname','arial')
tmp=pwd;
if isunix
    tmp2=max(find(tmp=='/'));
else
    tmp2=max(find(tmp=='\'));
    id_slash=find(tmp=='\');
    tmp2=id_slash(end-2);
end
tmp3=RemoveUnderScore(tmp((tmp2+1):end));
%tmp3=RemoveUnderScore(tmp((tmp2+1):end));
title_str{1}=['Site ID: ' tmp3 '     ' ineqal_str];
title_str{2}=' ';
title(title_str,'fontweight','bold')
axis off

if plot_spect_and_neural_ex
    axis on
    %    ax(1)=subplot(2,1,1);hold on;set(gca,'fontsize',12,'fontname','arial')

    if isunix
        ls *.cbin > batchfile
    else
        !dir /B *.cbin > batchfile
    end
    ls *.cbin
    fid=fopen('batchfile','r');
    true=0;
    while ~sum(true)
        fn=fgetl(fid);
        %disp('Skipping to 11th file');for zzz=1:10;        fn=fgetl(fid);end
        if (~ischar(fn))
            error(['Sequence "' search_for '" not found when looking for example .cbin file'])
            break;
        end
        if nargin==3        % plots example neural and song
            %            eval(sprintf('[dat,fs]=wavread(%s)',wav_file_name))
            disp(['Loading ' wav_file_name '.not.mat to look for correct label string x'])
            eval(sprintf('load %s.not.mat',wav_file_name))
            recfn=[fn(1:(end-4)) 'rec'];
            recdata=readrecf(recfn);        % .rec file contains timing of playbacks - FIRST file in BOS dir

            trigtime=recdata.ttimes(1)    % in msec

            onsets_neural = onsets+trigtime;    % these are different because neural and sound are from different files
            offsets_neural = offsets+trigtime;
            onsets_sound = onsets;
            offsets_sound = offsets;
            Fs_wav=Fs;
            Fs=32000;
        else
            disp(['Loading ' fn '.not.mat to look for correct label string'])
            eval(sprintf('load %s.not.mat',fn))
            onsets_neural = onsets;    % these are the SAME because neural and sound are from SAME files
            offsets_neural = offsets;
            onsets_sound = onsets;
            offsets_sound = offsets;
        end
        clear onsets offsets

        for x=1:(length(labels)-length(search_for)+1)
            compare=labels(x:(x+length(search_for)-1));
            true(x)=strcmp(search_for,compare);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        if timing_criterion    % change id_true to reflect some timing criterion
            id_true_search=find(true);
            for x=1:length(id_true_search)
                id=id_true_search(x):(id_true_search(x)+length(search_for)-1);
                next_interval_true(x,:)=(onsets_sound(id+1)-offsets_sound(id))/1000;
            end
            median_interval=median(next_interval_true(:,ID_interval));
            id_true_search_2=[];
            for x=1:length(id_true_search)
                if greater_than
                    if next_interval_true(x,ID_interval)>cutoff_interval
                        id_true_search_2=[id_true_search_2 id_true_search(x)];
                    end
                else
                    if next_interval_true(x,ID_interval)<cutoff_interval
                        id_true_search_2=[id_true_search_2 id_true_search(x)];
                    end
                end
            end
            disp([' Of ' num2str(length(id_true_search)) ' correct sequences, keeping ' num2str(length(id_true_search_2)) ' based on timing criterion'])
            true=zeros(1,length(onsets_sound));
            id=id_true_search_2;
            if ~isempty(id)
                if sum(id(1)==[16]);id=[];end
                %            if sum(id(1)==[16 43 49 150]);id=[];end
                %            if sum(id(1)==[142 152 43 14]);id=[];end
            end
            true(id)=1;
        else
            id=find(true);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        if sum(true)
            disp([ fn ' contains correct label string -  loading...'])
            eval(sprintf('load %s.neuralnot_CH%s.mat',fn,num2str(ch_with_data)))

            if nargin==3        % plots example neural and song
                disp([ 'Loading sound example from WAV file'])
                eval(sprintf('[song_sample,Fs_song_sample]=wavread(''%s'');',wav_file_name));
                th=1e-10;
            else
                disp([ 'Loading sound example from channel ' song_channel{1}])
                [song_sample, Fs_song_sample]=soundin('.', fn,song_channel{1});

                Fs_song_sample=32000;
                th=.1;
            end
            disp([ 'Loading neural example'])

            [neural_sample, Fs_neural_sample]=soundin('.', fn,['obs' num2str(ch_with_data)]);
            Fs_neural_sample=32000;
            if isfield(Data,'run_lowpass')
                if Data.run_lowpass
                    disp('Low-pass filtering example neural trace...')
                    neural_sample=bandpass(neural_sample,Fs_neural_sample,1,4000, 'hanningfir');
                end
            else
                disp('NOTE - lowpass filter field not yet added')
            end
            disp([ 'Finished loading'])
        end
    end
    disp('Starting spect')
    if use_preceding_interval
        clip_on_id=round(offsets_neural(id(1)-1)/1000*Fs);  % doesnt require an explicit vector of time
        clip_off_id=round(offsets_neural(id(1)+length(search_for)-1)/1000*Fs);
        clip_on_id_sound=round(offsets_sound(id(1)-1)/1000*Fs_song_sample);
        clip_off_id_sound=round(offsets_sound(id(1)+length(search_for)-1)/1000*Fs_song_sample);
    else
        clip_on_id=round(onsets_neural(id(1))/1000*Fs);  % doesnt require an explicit vector of time
        clip_off_id=round(offsets_neural(id(1)+length(search_for)-1)/1000*Fs);
        clip_on_id_sound=round(onsets_sound(id(1))/1000*Fs_song_sample);
        clip_off_id_sound=round(offsets_sound(id(1)+length(search_for)-1)/1000*Fs_song_sample);
    end

    neural_clip=neural_sample(clip_on_id:clip_off_id);

    %    L_time=length(clip_on_id_sound:clip_off_id_sound);
    L_time=length(clip_on_id:clip_off_id);
    %    neural_resize=neural_clip/max(abs(neural_clip));    % renorm to between roughly -1 and 1
    neural_resize=neural_clip;    % renorm to between roughly -1 and 1

    % invert neural waveform (if that's how spike time determination was
    % done)
    if Data.inverted_waveform
        neural_resize=-neural_resize;
    end

    [S1,F1,T1,P1] =spect_from_waveform(song_sample(clip_on_id_sound:clip_off_id_sound),Fs_song_sample,0,[0.8 8]); % will give a 512-point transform,
    %    [S1,F1,T1,P1] =spect_from_waveform(song_sample(clip_on_id:clip_off_id),Fs,0,[0.8 8]); % will give a 512-point transform,
    P1(P1<th)=th;
    displayspec_sam(T1,F1,P1,0,'yaxis');
    set(gca,'ylim',[500 10000])
    %set(gca,'Xlim',[0 offsets_true(1)-onsets_true(1)])
    tminmax_spect=[min(T1) max(T1)];
    set(gca,'Xlim',tminmax_spect);
    ymid=mean(get(gca,'ylim'));
    yrange=diff(get(gca,'ylim'));

    fclose(fid);

    spect_example_onsets=onsets_all_syls{1};
    spect_example_offsets=offsets_all_syls{1};

    ct=1;
    if use_preceding_interval
        spect_example_vec(ct)=spect_example_offsets(1);ct=ct+1;
        for x=1:length(spect_example_onsets);
            spect_example_vec(ct)=spect_example_onsets(x);ct=ct+1;  % canonical vec is on(1) off(1) on(2) off(2) etc
            spect_example_vec(ct)=spect_example_offsets(x+1);ct=ct+1;
        end
        spect_example_vec=spect_example_vec-spect_example_vec(1);   % set first time to zero
    else
        for x=1:length(spect_example_onsets);
            spect_example_vec(ct)=spect_example_onsets(x);ct=ct+1;  % canonical vec is on(1) off(1) on(2) off(2) etc
            spect_example_vec(ct)=spect_example_offsets(x);ct=ct+1;
        end
        spect_example_vec=spect_example_vec-spect_example_vec(1);   % set first time to zero
    end
    yl=get(gca,'ylim');
    for y=1:length(spect_example_vec)
        if rem(y,2)
            if y<(length(spect_example_vec))

                text(mean([spect_example_vec(y+use_preceding_interval) spect_example_vec(y+1+use_preceding_interval)]),yl(2)*1.1,search_for(round(y/2)),'horizontalalignment','center','fontsize',18,'fontweight','bold')
                %                text(mean([spect_example_vec(y) spect_example_vec(y+1)]),yl(2)*1.1,search_for(round(y/2)),'horizontalalignment','center','fontsize',18,'fontweight','bold')
            end
        end
    end
    ylabel('Freq (Hz)')

    disp('Spect done')
    ax(2)=subplot(8,1,3:4);hold on;set(gca,'fontsize',12,'fontname','arial')
    %    ydat=neural_resize(1:L_time)*.4*yrange+.5*yrange;
    ydat=neural_resize(1:L_time);
    xdat=[(0:L_time-1)/Fs];
    plot(xdat,ydat,'k','linew',1);
    set(gca,'ylim',[min(ydat) max(ydat)])
    set(gca,'xlim',tminmax_spect)
    axis off
    linkaxes(ax(1:2),'x')
end
xlabel('Time (sec)')
disp('Starting raster t=0');
vect_x=[];
vect_y=[];

for x=1:length(id_true)
    clip_on_id=round(onsets_true(x)*Fs);
    clip_off_id=round(offsets_true(x)*Fs);
    if isempty(clip_on_id) | isempty(clip_off_id)
        disp('empty')
    end

    ax(2)=subplot(2,1,2);hold on;set(gca,'fontsize',12,'fontname','arial')
    if ~isempty(spiketimes_warped{x})
        vect_x=[vect_x [spiketimes_warped{x};spiketimes_warped{x}]];
        vect_y=[vect_y [x*ones(1,length(spiketimes_warped{x})); (x+1)*ones(1,length(spiketimes_warped{x}))]];
        st_quants=floor(spiketimes_warped{x}/tbin_for_hist)+1;
        quantbins=unique(st_quants);
        [b,c]=histc(st_quants,quantbins);
        MAT(x,quantbins)=b;
    end

    if plot_spect_and_neural_ex & x==1
        t_points=spiketimes_warped{x};
        xxx=gca;
        subplot(8,1,3:4);
        yl=get(gca,'ylim');
        xl=get(gca,'xlim');
        plot(first_spiketrain_unwarped,ydat(round(first_spiketrain_unwarped*32000)+1),'ko','markerfacecolor',spike_color,'markersize',5)
        %        plot(first_spiketrain_unwarped,ydat(round(first_spiketrain_unwarped*32000)+1),'ko','markerfacecolor','r','markersize',5)

        set(gcf,'currentaxes',xxx)
    end

end

set(gca,'xlim',[0 max(canonical_vec)])
disp('Raster done')

rate_out=full(sum(MAT,1)/(tbin_for_hist *length(id_true)));
std_out=full(std(MAT)/tbin_for_hist);
cv_out=std_out./rate_out;

m_rate=round(max(rate_out));

[AX,H1,H2] =  plotyy((0:(size(MAT,2)-1))*tbin_for_hist,rate_out,vect_x,vect_y);  % poor man's raster plot
set(AX(1),'xlim',[0 max(canonical_vec)],'ycolor','k')
%set(H1,'color','b','linew',2)
set(H1,'color',spike_color,'linew',2)
set(gcf,'currentaxes',AX(1));hold on
ylabel(AX(1),'Firing rate (Hz)')

if sort_by_pitch
    set(H1,'visible','off')
    set(gcf,'currentaxes',AX(2));hold on
    if mark_intervals_and_pitchquant
        for x=1:length(id_true)
            plot(spike_intervals_warped{x},[x x],'ro','markerfacecolor','r','markersize',2)
        end
    end
end

set(AX(2),'xlim',[0 max(canonical_vec)],'ycolor','k','ytick',[1 length(id_true)+1],'ylim',[1 length(id_true)+1],'fontsize',12)
%set(AX(2),'xlim',[0 max(canonical_vec)],'ycolor','k','ytick',[1 max(max(vect_y))],'ylim',[1 max(max(vect_y))],'fontsize',12)
if sort_by_pitch
    ylabel(AX(2),'# trials, sorted by pitch')
else
    ylabel(AX(2),'# trials')
end
clear tmp
tmp{1}='  ';
tmp{2}=' Time (sec)';
xlabel(tmp)
set(AX(1),'xtick',[])
set(H2,'color','k')%
%linkaxes([ax AX],'x');
linkaxes([AX],'x');

%return

if 0
    linkaxes(ax,'x');
    rate_out=full(sum(MAT)/(tbin_for_hist *length(id_true)));

    m_rate=round(max(rate_out));
    rate_out=rate_out/max(rate_out)*(length(id_true))+1;
    set(gca,'ylim',[1 length(id_true)+1])
    set(gca,'xlim',tminmax_spect)
    set(gca,'ytick',yl)
    ylab{1}=num2str(0);
    ylab{2}=num2str(m_rate);
    plot((0:(size(MAT,2)-1))*tbin_for_hist,rate_out,'w','linew',4)
    plotyy(vect_x,vect_y,(0:(size(MAT,2)-1))*tbin_for_hist,rate_out);  % poor man's raster plot
    set(gca,'Xlim',tminmax_spect)
    tmp{1}=' ';
    tmp{2}=' Time (sec)';
    xlabel(tmp)
    ylabel('Firing rate (Hz)')
end



if lines_connecting_songs % draw lines connecting songs from same file
    set(gcf,'currentaxes',AX(2));hold on
    xl=get(gca,'xlim');
    yl=get(gca,'ylim');
    [a,b]=unique(sequence_from_file_number);
    %vline_lims=cumsum([1 b]);
    vline_lims=[0 b];
    for x=2:length(vline_lims)
        plot([xl(2) xl(2)],[vline_lims(x-1)+1.25 vline_lims(x)+.75],'r','linew',5)
        tmp=fname_arr{a(x-1)};
        id_tmp=max(find(tmp=='_'))+1;
        use_txt=RemoveUnderScore(tmp(id_tmp:end));
        text(xl(2),mean([vline_lims(x-1)+1.25 vline_lims(x)+.75]),[' ' use_txt],'horizontalalignment','left','color','r')
    end
else
    if use_preceding_interval
        canonical_vec=canonical_vec(2:end);
    end
    plot_pink_boxes(canonical_vec,search_for)
end

function plot_pink_boxes(canonical_vec,search_for)

xl=get(gca,'xlim');
yl=get(gca,'ylim');
for y=1:length(canonical_vec)
    if rem(y,2)
        if y<(length(canonical_vec))
            text(mean([canonical_vec(y) canonical_vec(y+1)]),yl(2)*1.1,search_for(round(y/2)),'horizontalalignment','center','fontsize',18,'fontweight','bold')
        end
        if 1    % draws pink boxes for syllables - problem is that these only render right using the openGL renderer.
            % but if pink boxes included, cant get to print as anything but
            % image.  only painter's uses vector graphics
            x_pts=[canonical_vec(y) canonical_vec(y)  canonical_vec(y+1) canonical_vec(y+1)];
            y_pts=[yl fliplr(yl)];
            %        a=patch(x_pts,y_pts,'w');
            a=fill(x_pts,y_pts,'w');
            set(a,'facecolor',[1 .9 .9])
            set(a,'edgecolor',[1 .9 .9])
            set(gcf,'renderer','opengl')    % screen
            %            set(gcf,'renderer','painters') % printed fig
            %        set(a,'facecolor',[.9 .9 1])
            %        set(a,'facecolor',[.8 .8 .8])
            %            set(a,'edgealpha',0)
        end
    end
end

function  out=compute_summed_amp_window(smoothed,birdname,labels);

if strcmp(birdname,'B53O71') | strcmp(birdname,'O55Pu53') | strcmp(birdname,'G26G23') ...
        | strcmp(birdname,'G34G54')    | strcmp(birdname,'G45G46') | strcmp(birdname,'Pu55Pu22')
    [f_cutoff,t_assay,spect_params]=syllable_params_by_bird_MEL(birdname,labels);
else
    [f_cutoff,t_assay,spect_params]=syllable_params_by_bird(birdname,labels);
end
window_half_width=8;   % in msec, the width of the window on EITHER side of time_quant
%window_half_width=12.5   % in msec, the width of the window on EITHER side of time_quant
time=[0:length(smoothed)-1]/32000;
if strcmp(f_cutoff,'undefined')
    out=0;
else
    if t_assay>1    % if t_assay>1, then it is the percentage (1-100) of the syl at which to quantify
        time_quant=max(time)*t_assay*.01;
    else
        time_quant=t_assay;
    end
    window=time_quant+[-window_half_width window_half_width]*.001;% in seconds
    window_id=round(window*32000);% id of time at start and end of window
    if window_id(1)<1;window_id(1)=1;end
    if window_id(2)>length(smoothed);window_id(2)=length(smoothed);end
    out=mean(smoothed(window_id(1):window_id(2)));
    if isnan(out);
        labels
        t_assay
        window_id,
%        error('NaN');
out=0;disp('NaN for amplitude, setting to zero')
    end
    if length(out)>1;error('wer');end
    if 0%labels=='n'
        figure(100);clf;hold on
        plot(time,log10(smoothed))
        yl=get(gca,'ylim');
        plot([time_quant time_quant],yl,'k:')
        plot([window(1) window(1)],yl,'g:')
        plot([window(2) window(2)],yl,'r:')
        title(num2str(out))
        labels
        pause(3)
    end
end
