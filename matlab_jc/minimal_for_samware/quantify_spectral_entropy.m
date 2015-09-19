%
% copied from quantify_pitch.m
%
% meant to be used on phys only data.  adds computed peak pitch to .not.mat
%

function foo(syl_to_quant,channel_with_song)

if strcmp(channel_with_song,'obs0r')
    ch_with_song='obs0r';disp(' ');disp('Loading song from channel 0r');disp(' ')
elseif strcmp(channel_with_song,'obs0')
    ch_with_song='obs0';
end

if isunix
    !rm -f batchfile
    !ls -1 *.not.mat > batchfile
else
    !dir /B **.not.mat > batchfile
end
if strfind(pwd,'pu44w52');bird_name='pu44w52';
elseif strfind(pwd,'pu24w39');bird_name='pu24w39';
elseif strfind(pwd,'g26g23');bird_name='g26g23';
elseif strfind(pwd,'pu24w39');bird_name='pu24w39';
elseif strfind(pwd,'r75g59');    bird_name='r75g59';
elseif strfind(pwd,'pu5b3');    bird_name='pu5b3';
elseif strfind(pwd,'w48o28');    bird_name='w48o28';
elseif strfind(pwd,'pu54b57');    bird_name='pu54b57';
elseif strfind(pwd,'pu53b58');    bird_name='pu53b58';
elseif strfind(pwd,'pu28b75');    bird_name='pu28b75';
elseif strfind(pwd,'g18g8');    bird_name='g18g8';
    % NCM implant
elseif strfind(pwd,'ZFpi90w37');    bird_name='ZFpi90w37';
    % NCM implant
elseif strfind(pwd,'g44o23');    bird_name='g44o23';
elseif strfind(pwd,'pu26y2');    bird_name='pu26y2';
elseif strfind(pwd,'r12r11');    bird_name='r12r11';
    divide_syl_target{1}='c';          % replaces target syl with _replacement.  use this if want to quantify pitch at two place
    divide_syl_replacement{1}='mn'; % in one syl
elseif strfind(pwd,'G26G23');    bird_name='G26G23';
elseif strfind(pwd,'Pu55Pu22');    bird_name='Pu55Pu22';
elseif strfind(pwd,'W90O73');    bird_name='W90O73';
elseif strfind(pwd,'O55Pu53');    bird_name='O55Pu53';
elseif strfind(pwd,'O55Pu53');    bird_name='O55Pu53';
elseif strfind(pwd,'W32Pi51');    bird_name='W32Pi51';
elseif strfind(pwd,'W15W94');    bird_name='W15W94';
elseif strfind(pwd,'O14O15');    bird_name='O14O15';
elseif strfind(pwd,'G45G46');    bird_name='G45G46';
elseif strfind(pwd,'Pk35G27');    bird_name='Pk35G27';
elseif strfind(pwd,'B53O71');    bird_name='B53O71';
elseif strfind(pwd,'W96Pi45');    bird_name='W96Pi45';
elseif strfind(pwd,'g38o18');    bird_name='g38o18';
elseif strfind(pwd,'p85g54');    bird_name='p85g54';
elseif strfind(pwd,'b39b14');    bird_name='b39b14'; % mimi lman ZF
elseif strfind(pwd,'p49w84');    bird_name='p49w84'; % mimi lman ZF
elseif strfind(pwd,'o99b51');    bird_name='o99b51'; % mel lman
elseif strfind(pwd,'b9r63');    bird_name='b9r63'; % mimi lman ZF
elseif strfind(pwd,'blueblue');    bird_name='blueblue'; % mimi lman ZF
elseif strfind(pwd,'jman3');    bird_name='jman3'; % mimi lman ZF
elseif strfind(pwd,'masa');    bird_name='masa'; % mimi lman ZF
elseif strfind(pwd,'pur98w3');    bird_name='pur98w3'; % mimi lman ZF
elseif strfind(pwd,'pu77bk41');    bird_name='pu77bk41'; % mimi lman ZF
elseif strfind(pwd,'o85pu54');    bird_name='o85pu54'; % mimi lman ZF
elseif strfind(pwd,'g91pu54');    bird_name='g91pu54';
elseif strfind(pwd,'bl82bl81');    bird_name='bl82bl81';
elseif strfind(pwd,'r93bl81');    bird_name='r93bl81';
end


disp(['Bird name = ' bird_name])
fid_extra=fopen('batchfile','r');
fn_extra=fgetl(fid_extra);
uid=find(fn_extra=='_');
if isempty(bird_name)   % if wav files
    uid=find(fn_extra=='.');bird_name=fn_extra(1:uid-1)
end
fclose(fid_extra);
disp('calling syllable_params_by_bird')

if  strcmp(bird_name,'G26G23') |  strcmp(bird_name,'Pu55Pu22') |  strcmp(bird_name,'W90O73') |  strcmp(bird_name,'O55Pu53') |  strcmp(bird_name,'W32Pi51') |  strcmp(bird_name,'W15W94') |  strcmp(bird_name,'O14O15') |  strcmp(bird_name,'G45G46')|  strcmp(bird_name,'Pk35G27')|  strcmp(bird_name,'B53O71')|  strcmp(bird_name,'W96Pi45')
    [f_cutoff,t_assay,spect_params]=syllable_params_by_bird_MEL(bird_name,syl_to_quant);
    syl_to_quant
else
    [f_cutoff,t_assay,spect_params]=syllable_params_by_bird(bird_name,syl_to_quant);
end
if t_assay>1    % if t_assay>1, then it is the percentage (1-100) of the syl at which to quantify
    use_pct=1;t_pct=t_assay*.01;
    disp(['Assaying at ' num2str(t_assay) '% of total syllable duration with freq cutoff at ' num2str(f_cutoff)])
else
    use_pct=0;    disp(['Assaying at t=' num2str(t_assay) ' with freq cutoff at ' num2str(f_cutoff)])
end
disp(' ');disp(['Window duration ' num2str(spect_params(2)) ' msec, overlap = ' num2str(spect_params(1)) ]);disp(' ')

fid=fopen('batchfile','r');
save_T1=[];
ct_nfiles=1;

while (1)

    fn=fgetl(fid);
    if ct_nfiles==1;
        uid=find(fn=='_');
        %        bird_name=fn(1:uid-1);
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
    clear spectral_entropy
    eval(sprintf('load %s',fn))
    ct=1;

    eval(sprintf('field_summary=load(''%s'');',fn))
    % BELOW NEEDS TO COME AFTER LOADING OR IT WILL BE OVERWRITTEN
    % trim freq and power to range of song power
    %    power_range_for_spect_entropy=[500 10000];


    fieldnames_arr=fields(field_summary);
    fieldnames=[];for x=1:length(fieldnames_arr);fieldnames=[fieldnames ' ' fieldnames_arr{x}];end
    if strcmp(cbin_fn(end-3:end),'cbin')
        [rawsong, Fs]=evsoundin('.', cbin_fn,ch_with_song);
        %        [rawsong, Fs]=evsoundin('.', cbin_fn,ch_with_song);
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

    Ls=labels;

    % get onsets and offsets into sec, not msec
    onsets=onsets/1000;offsets=offsets/1000;

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
        if length(T1)>length(save_T1);save_T1=T1;end

        %         OLD WAY
        % %         % trim freq and power to range of song power
        % %         f_cut_id=find(F1>power_range_for_spect_entropy(1) & F1<power_range_for_spect_entropy(2));
        % %         F1=F1(f_cut_id);
        % %         P1=P1(f_cut_id,:);
        % %         t_id=find(abs((T1-t_assay))==min(abs(T1-t_assay)));
        % %         spect_slice=P1(:,t_id);

        % NEW WAY
        [f_cutoff,t_assay,spect_params]=syllable_params_by_bird(bird_name,syl_to_quant);
        if strcmp(spect_params,'undefined')
            [f_cutoff,t_assay,spect_params]=syllable_params_by_bird_MEL(bird_name,syl_to_quant);
            disp('Using syllable_params_by_bird_MEL.m')
        end
        f_cut_id_tmp=find(F1>f_cutoff(1) & F1<f_cutoff(2));
        F1_tmp=F1(f_cut_id_tmp);
        P1_tmp=P1(f_cut_id_tmp,:);

        t_id=find(abs((T1-t_assay))==min(abs(T1-t_assay)));
        spect_slice_tmp=P1_tmp(:,t_id(1));
        
        max_p_id_tmp=find(spect_slice_tmp==max(spect_slice_tmp));
        peak_tmp=F1_tmp(max_p_id_tmp);

        % %         f_cut_id=find(F1>peak_tmp*2^-.5 & F1<peak_tmp*2^.5);    % one octave of freqs, centered on peak
        % %        disp([ F1(f_cut_id(1))        F1(f_cut_id(end))])

        harm_number=1%get_harmonic_number(syl_to_quant);
        f_1=peak_tmp/harm_number;
        %         cut=f_1*2^-.5;
        %         range=cut*[2^(harm_number-1) 2^harm_number];
        %         disp(['peak_tmp = ' num2str(peak_tmp)])
        %         disp(['f_1 = ' num2str(f_1)])
        %         disp(['range = ' num2str(range)])

        % range equal to half of f_1, centered (numerically) on peak
        %alt_range=peak_tmp+[-f_1/4 +f_1/4];
        % range equal to  f_1, centered (numerically) on peak
        alt_range=peak_tmp+[-f_1/2 +f_1/2];
        disp(['alt_range, peak, f1 = ' num2str([alt_range peak_tmp f_1])])
        range=alt_range;%disp('USING ALT RANGE')

        %         alt_cut=cut/2;
        %         disp(['cut = ' num2str(cut)])
        %         disp(['alt_cut = ' num2str(alt_cut)])
        %         alt_range=alt_cut*([2^(harm_number) 2^(harm_number+1)]);
        %         disp(['alt_range = ' num2str(alt_range)])

        f_cut_id=find(F1>range(1) & F1<range(2));    % one octave of freqs, centered on peak



        F1=F1(f_cut_id);
        P1=P1(f_cut_id,:);


        t_id=find(abs((T1-t_assay))==min(abs(T1-t_assay)));
        spect_slice=P1(:,t_id);

        technique=3;

        if technique==1
            disp('THRESHOLDING POWER')
            power_TH=3;        %        power_TH=300;
            spect_slice_TH=spect_slice(find(spect_slice>power_TH));
            norm_spect_TH=spect_slice_TH/sum(spect_slice_TH);
            spectral_entropy_tmp(ct)=-sum(norm_spect_TH.*log10(norm_spect_TH));

            % for plotting
            norm_spect=spect_slice/sum(spect_slice);
        elseif technique==2
            disp('SQRTing POWER')
            norm_spect=spect_slice.^.5/sum(spect_slice.^.5);
            spectral_entropy_tmp(ct)=-sum(norm_spect.*log10(norm_spect));
        elseif technique==3
            norm_spect=spect_slice/sum(spect_slice);
            spectral_entropy_tmp(ct)=-sum(norm_spect.*log10(norm_spect));
        end

        make_plots=0;
        if make_plots
            ofst=1600;
            % pu24w39
            range_1=[.8 .9];  %syl d
            range_2=[.4 .5];

            %pu26y2
            %             range_1=[.8 1.2]; % syl b, alt range
            %             range_2=[.6 .7];
            %             range_1=[.8 1.1];% syl b
            %             range_2=[.5 .6];
            %             range_1=[.8 .9];% syl c
            %             range_2=[.4 .5];
            %             range_1=[.8 1]; % syl d
            %             range_2=[.4 .6];

            if spectral_entropy_tmp(ct)>range_1(1) & spectral_entropy_tmp(ct)<range_1(2)
                figure(400+ofst);ax(1)=subplot(1,2,1);hold on;plot(F1,spect_slice,'r');title(['HIGH Entropy = ' num2str(spectral_entropy_tmp(ct))])
                set(gcf,'paperposition',[0.2500    2.5000    8.0000    2.8656])
                figure(401+ofst);subplot(2,1,1);hold on;plot(F1,log10(spect_slice),'r'),title('Red is high entropy, blue is low')
                figure(401+ofst);subplot(2,1,2);hold on;plot(F1,-(norm_spect.*log10(norm_spect)),'r'),title('Red is high entropy, blue is low')
                figure(402+ofst);bx(1)=subplot(2,1,1);hold on;spect_from_waveform(syl_wav,Fs,1,spect_params);title('high entropy')
                set(gca,'ylim',[500 10000]); yl=get(gca,'ylim');plot(t_assay*[1 1],yl,'k')
                set(gcf,'paperposition',[0.2500    2.5000    1.8375    6.0000])
            elseif spectral_entropy_tmp(ct)>range_2(1)  & spectral_entropy_tmp(ct)<range_2(2)
                figure(400+ofst);ax(2)=subplot(1,2,2);hold on;plot(F1,spect_slice);title(['LOW Entropy = ' num2str(spectral_entropy_tmp(ct))])
                set(gcf,'paperposition',[0.2500    2.5000    8.0000    2.8656])
                figure(401+ofst);subplot(2,1,1);;hold on;plot(F1,log10(spect_slice),'b')
                figure(401+ofst);subplot(2,1,2);hold on;plot(F1,-(norm_spect.*log10(norm_spect)),'b'),
                set(gcf,'paperposition',[2.4281    2.9250    3.1125    6.0000])
                figure(402+ofst);bx(2)=subplot(2,1,2);hold on;spect_from_waveform(syl_wav,Fs,1,spect_params);title('low entropy')
                set(gca,'ylim',[500 10000]);yl=get(gca,'ylim');plot(t_assay*[1 1],yl,'k')
                set(gcf,'paperposition',[0.2500    2.5000    1.8375    6.0000])
            end
        end
        fname_arr{ct}=cbin_fn;
        ida=max(find(cbin_fn=='_'))+1;    % last character '_' +1
        idb=min(find(cbin_fn=='.'))-1;    % first character '.' -1


        ct=ct+1;
    end

    if ~exist('spectral_entropy')
        spectral_entropy=zeros(size(labels));
        disp('No previous spectral_entropy detected')
    elseif length(spectral_entropy)<length(labels)
        spectral_entropy=zeros(size(labels));
        disp('spectral_entropy not long enough - replaced by zeros')
    end

    ctX=1;
    for x=find(labels==syl_to_quant)
        spectral_entropy(x)=spectral_entropy_tmp(ctX);
        ctX=ctX+1;
    end
    onsets=onsets*1000;offsets=offsets*1000;
    if exist('range')==2;range=[0 0];end
    eval(sprintf('save %s %s',fn,[fieldnames ' range  channel_with_song spectral_entropy']))
end
if make_plots
    linkaxes(ax,'xy')
    linkaxes(bx,'xy')
end
