%
% copied from time_vs_f1peak.m
%
% meant to be used on phys only data.  adds computed peak pitch to .not.mat
%

function foo(syl_to_quant,channel_with_song)

if strcmp(channel_with_song,'obs0r')
ch_with_song='obs0r';disp(' ');disp('Loading song from channel 0r');disp(' ')
elseif strcmp(channel_with_song,'obs0')
ch_with_song='obs0';
end
!rm -f batchfile
if isunix
    !ls -1 *.not.mat > batchfile
else
    !dir /B **.not.mat > batchfile

    %    !dir /B *1728*.not.mat > batchfile
end
if strfind(pwd,'pu44w52')
    bird_name='pu44w52';
elseif strfind(pwd,'g26g23')
    bird_name='g26g23';
elseif strfind(pwd,'pu24w39')
    bird_name='pu24w39';
elseif strfind(pwd,'r75g59')
    bird_name='r75g59';
elseif strfind(pwd,'pu5b3')
    bird_name='pu5b3';
elseif strfind(pwd,'w48o28')
    bird_name='w48o28';
elseif strfind(pwd,'pu54b57')
    bird_name='pu54b57';
elseif strfind(pwd,'pu53b58')
    bird_name='pu53b58';
elseif strfind(pwd,'pu28b75')
    bird_name='pu28b75';
elseif strfind(pwd,'g18g8')
    bird_name='g18g8';
elseif strfind(pwd,'ZFpi90w37')% NCM implant
    bird_name='ZFpi90w37';
elseif strfind(pwd,'g44o23')% NCM implant
    bird_name='g44o23';
elseif strfind(pwd,'pu26y2')
    bird_name='pu26y2';
elseif strfind(pwd,'r12r11')
    bird_name='r12r11';
    divide_syl_target{1}='c';          % replaces target syl with _replacement.  use this if want to quantify pitch at two place 
    divide_syl_replacement{1}='mn'; % in one syl
elseif strfind(pwd,'pu77bk41')
    bird_name='pu77bk41';
elseif strfind(pwd,'o85pu54')
    bird_name='o85pu54';
elseif strfind(pwd,'g91pu54')
    bird_name='g91pu54';
elseif strfind(pwd,'bl82bl81')
    bird_name='bl82bl81';
elseif strfind(pwd,'r93bl81')
    bird_name='r93bl81';
end


bird_name
fid_extra=fopen('batchfile','r');
fn_extra=fgetl(fid_extra);
uid=find(fn_extra=='_');
%bird_name=fn_extra(1:uid-1)
if isempty(bird_name)   % if wav files
    uid=find(fn_extra=='.');
    bird_name=fn_extra(1:uid-1);
end
fclose(fid_extra);
disp('calling syllable_params_by_bird')

[f_cutoff,t_assay,spect_params]=syllable_params_by_bird(bird_name,syl_to_quant)

if t_assay>1    % if t_assay>1, then it is the percentage (1-100) of the syl at which to quantify
    use_pct=1;
    t_pct=t_assay*.01;
    disp(['Assaying at ' num2str(t_assay) '% of total syllable duration with freq cutoff at ' num2str(f_cutoff)])
else
    use_pct=0;
    disp(['Assaying at t=' num2str(t_assay) ' with freq cutoff at ' num2str(f_cutoff)])
end

disp(' ')
disp(['Window duration ' num2str(spect_params(2)) ' msec, overlap = ' num2str(spect_params(1)) ])
disp(' ')



fid=fopen('batchfile','r');
save_T1=[];
ct_nfiles=1;
%ct=1;
while (1) %& ct<3
    fn=fgetl(fid);
    if ct_nfiles==1;
        uid=find(fn=='_');
        bird_name=fn(1:uid-1);
        ct_nfiles=ct_nfiles+1;
    end
    cbin_fn=fn(1:end-8);
    rec_fn=[cbin_fn(1:end-5) '.rec'];
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    disp([fn ' quantifying syllable ' syl_to_quant]);
    clear peak_labelvec peak_pinterp_labelvec amplitude_at_pitchquant sum_power peak peak_pinterp t_out
    eval(sprintf('load %s',fn))
    ct=1;

    
    eval(sprintf('field_summary=load(''%s'');',fn))
    fieldnames_arr=fields(field_summary);
    fieldnames=[];for x=1:length(fieldnames_arr);fieldnames=[fieldnames ' ' fieldnames_arr{x}];end
    if strcmp(cbin_fn(end-3:end),'cbin')
        [rawsong, Fs]=evsoundin('.', cbin_fn,ch_with_song);
        is_cbin=1;
    else    % if file is .wav
        [rawsong,Fs,NBITS]=wavread(cbin_fn);
        is_cbin=0;
        rawsong=resample(rawsong,32000,Fs);
    end
    Fs=32000;
    
    if exist('divide_syl_target')
        labels_old=labels; labels=[];
        onsets_old=onsets;onsets=[];
        offsets_old=offsets;offsets=[];
        for x=1:length(labels_old)
            for y=1:length(divide_syl_target)
                if ~sum(labels_old(x)==divide_syl_target{y})
                    labels=[labels labels_old(x)];
                    onsets=[onsets; onsets_old(x)'];
                    offsets=[offsets; offsets_old(x)'];
                else
                    labels=[labels divide_syl_replacement{y}];                    
                    onsets=[onsets; onsets_old(x)*ones(length(divide_syl_replacement{y}),1)];                    
                    offsets=[offsets; offsets_old(x)*ones(length(divide_syl_replacement{y}),1)];                    
                    end
            end
        end
    end
%     labels_old
% labels
% return
    Ls=labels;

    % get onsets and offsets into sec, not msec
    onsets=onsets/1000;
    offsets=offsets/1000;

    t_samples=1:length(rawsong);
    t=t_samples/Fs;

    id=find(Ls==syl_to_quant);
    disp(['Found ' num2str(length(id)) ' examples'])
    for x=1:length(id)

        on=onsets(id(x));
        off=offsets(id(x));
        if use_pct
            t_assay=t_pct*(off-on);
        end
        if (off-on)<spect_params(2)/1000    %if syl is too short
            off=on+spect_params(2)/1000;    % extend end of syl
            disp('Extending offset')
        end
        on_id=find(abs((t-on))==min(abs(t-on)));
        off_id=find(abs((t-off))==min(abs(t-off)));
        syl_wav=rawsong(on_id:off_id);

        [S1,F1,T1,P1] =spect_from_waveform(syl_wav,Fs,0,spect_params);
        %disp(['Time to run spect - ' num2str(toc)])
        P1_uncut=P1;
        F1_uncut=F1;
        [r,c]=size(P1);
        whole_P1(1:r,1:c,ct)=P1;
        whole_F1{ct}=F1;
        if length(T1)>length(save_T1)
            save_T1=T1;
        end
        % TO BE REPLACED BY FUN
        f_cut_id=find(F1>f_cutoff(1) & F1<f_cutoff(2));
        F1=F1(f_cut_id);
        P1=P1(f_cut_id,:);


        t_id=find(abs((T1-t_assay))==min(abs(T1-t_assay)));
        spect_slice=P1(:,t_id);
        slice_save{ct,1}=spect_slice;

        %figure(2);clf;plot(F1_uncut,P1_uncut(:,t_id));return
        % a measure of amplitude - sum power spectrum
        sum_power(ct,1)=sum(P1_uncut(:,t_id));
%        plot(whole_F1{ct},P1_uncut(:,t_id)),return

        F1_save{ct}=F1;
        max_p_id=find(spect_slice==max(spect_slice));
        peak(ct,1)=F1(max_p_id);
        

        if sum(max_p_id==[1 length(spect_slice)])
        if max_p_id==1
            disp('Peak at LOWER extreme of range')
        elseif max_p_id==length(spect_slice)
            disp('Peak at UPPER extreme of range')
        end
            peak_pinterp(ct,1)=peak(ct,1);
        else
            x_pinterp=F1(max_p_id-1:max_p_id+1);
            y_pinterp=P1(max_p_id-1:max_p_id+1,t_id);
            [peak_pinterp(ct,1),tmp] = pinterp(x_pinterp, y_pinterp);
            if peak_pinterp(ct,1)==0
                disp(['zero here     ' fn])
            end
        end



        fname_arr{ct}=cbin_fn;
        ida=max(find(cbin_fn=='_'))+1;    % last character '_' +1
        idb=min(find(cbin_fn=='.'))-1;    % first character '.' -1


        % calculate time from .REC header
        if is_cbin
            recdata=readrecf(rec_fn);
            t_header_str=recdata.header{1};
            space_id=find(t_header_str==' ');
            t_str=t_header_str((space_id(end)+1):end);
            hr=str2num(t_str(1:2));
            minute=str2num(t_str(4:5));
            second=str2num(t_str(7:8));
            t_out(ct)=hr+minute/60+second/3600+on/3600;
        else
            t_out(ct)=ct;
        end

        idc=min(find(cbin_fn=='_'))+1;    % first character '_' +1

        ct=ct+1;
        % get a string for the date and syllable
        ul_id=find(cbin_fn=='_');
        if numel(ul_id) > 1
%        if ~isempty(ul_id) & length(ul_id)<=1
            if is_cbin
                idd=ul_id(end-1);
                date_str=cbin_fn(idd+1:(ida-2));
%                date_str=[date_str(3:4) '_' date_str(1:2) '_20' date_str(5:6) '_syl_' syl_to_quant];
            else
                p_id=find(cbin_fn=='.');
                pdd=p_id(1);
                date_str=cbin_fn(pdd+1:pdd+8);
                date_str=[date_str(5:6) '_' date_str(7:8) '_' date_str(1:4) '_syl_' syl_to_quant '_WAV'];
            end
        else
            date_str='NODATE';
            if ct==1
                disp('Date string not extracted from filename')
            end
        end
    end

    if ~exist('amplitude_at_pitchquant')
        amplitude_at_pitchquant=zeros(size(labels));
        disp('No previous amplitude_at_pitchquant detected')
    elseif length(amplitude_at_pitchquant)<length(labels)
        amplitude_at_pitchquant=zeros(size(labels));
        disp('amplitude_at_pitchquant not long enough - replaced by zeros')
    end
    if ~exist('peak_labelvec')
        peak_labelvec=zeros(size(labels));
        peak_pinterp_labelvec=zeros(size(labels));
        disp('No previous peaklabels detected')
    else
        disp('Overwriting/adding to previously computed peaks')
    end
    ctX=1;
    for x=find(labels==syl_to_quant)
        amplitude_at_pitchquant(x)=sum_power(ctX);
        peak_labelvec(x)=peak(ctX);
        peak_pinterp_labelvec(x)=peak_pinterp(ctX);
        ctX=ctX+1;
    end
    onsets=onsets*1000;
    offsets=offsets*1000;
    eval(sprintf('save %s %s',fn,[fieldnames ' peak_labelvec peak_pinterp_labelvec amplitude_at_pitchquant channel_with_song']))
end
