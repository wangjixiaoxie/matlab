function t_peaks=foo(file,ch,analyze);
if nargin==2
    analyze=0;
    batchmode=0;
end
if nargin==3
    batchmode=0;
end
if nargin==1
    batchmode=1;
end

if ~batchmode
    analyze_func(file,ch,analyze)
else
    eval(sprintf('load %s',file));%batchmodefile should be a .neuralnote.mat file
    disp(' ')
    disp(['Running in batchmode.'])
    disp(['Using setting loaded from ' file])
    disp(['Data.channels = ' num2str(Data.channels)])
    disp(['Data.inverted_waveform = ' num2str(Data.inverted_waveform)])
    disp(['Data.TH = ' num2str(Data.TH)])

    disp(' ')
    if isunix
        !ls *.cbin > batchfoo
    else
        !dir /B *.cbin > batchfoo
    end
    fid=fopen('batchfoo','r');
    while 1
        fn=fgetl(fid);
        if (~ischar(fn))
            break;
        end
        %        file=fn(1:end-14);% cuts off '.neuralnot.mat'
        analyze_func(fn,Data.channels,1:4,Data)  % Data provides parameters for doing analysis
    end
    fclose(fid);
end

function analyze_func(file,ch,analyze,Data)
if nargin<4
    batchmode=0;
else
    batchmode=1;
    invert=Data.inverted_waveform;
    low_filt=Data.lowpass_filter;
    TH=Data.TH;
    %ch_select=Data.channels;

    ch_tmp=Data.channels;
    ch_select=intersect(ch_tmp,ch);
    for x=1:length(ch)
        z(x)=sum(ch(x)==ch_tmp);
    end
    ch_select=find(z);

    clear Data
end
TH_default=400;
%TH_default=1e4;

if ~batchmode
    inv_str=input('Invert waveforms? [y/n], or a vector of channels to invert - [x x x] ','s');
    if strcmp(inv_str,'y')
        invert=ones(1,length(ch));
    elseif strcmp(inv_str,'n')
        invert=zeros(1,length(ch));
    else
        invert=zeros(1,length(ch));
        inv_id=str2num(inv_str);
        for x=1:length(ch)
            if sum(inv_id==ch(x))
                invert(x)=1;
            end
        end
    end
end

for x=1:length(ch)
    if ch(x)==1;
        ch_str_array{x}='obs1';
    elseif ch(x)==2;
        ch_str_array{x}='obs2';
    elseif ch(x)==3;
        ch_str_array{x}='obs3';
    elseif ch(x)==4;
        ch_str_array{x}='obs4';
    elseif ch(x)==10;
        ch_str_array{x}='obs1r';
    elseif ch(x)==20;
        ch_str_array{x}='obs2r';
    elseif ch(x)==30;
        ch_str_array{x}='obs3r';
    elseif ch(x)==40;
        ch_str_array{x}='obs4r';
    end
    disp(['loading ' file ', channel ' ch_str_array{x}])
    [neural{x}, Fs]=evsoundin('.', file,ch_str_array{x});
%    [neural{x}, Fs]=soundin('.', file,ch_str_array{x});
    Fs=32000;
    if sum(invert)
        if invert(x)
            neural{x}=-1*neural{x};
            disp(['inverting channel ' ch_str_array{x}])
        end
    end
end


if batchmode & sum(low_filt)
    for x=1:length(ch)
        if   low_filt(x);
            disp(['Low-pass filtering channel ' num2str(ch(x)) ' to remove movement artifact']);
            neural{x}=bandpass(neural{x},Fs,1,4000, 'hanningfir');
        end
    end
end

if ch(1)<6
    [rawsong, tmp]=evsoundin('.', file,'obs0');
else
    [rawsong, tmp]=evsoundin('.', file,'obs0r');
end

if Fs==3;Fs=32000;end   % I still dont understant this bug
tmax=length(neural{1})/Fs;
time=(1/Fs):(1/Fs):tmax;


if ~batchmode
    nfigs=1+length(ch);
    a=gcf;
    ax(1)=subplot(nfigs,1,1);

    [S1,F1,T1,P1] =spect_from_waveform(rawsong,Fs,0,[0.8 8]); % will give a 512-point transform,
    % recall that CBS has a FS=44100
    P1(P1<.001)=.001;
    displayspec_sam(T1,F1,P1,0,'yaxis')
    set(gca,'ylim',[500 10000])


    title([file ', sound channel obs0'])
    for x=1:length(ch)
        ax(x+1)=subplot(nfigs,1,x+1);hold on
        plot(time,neural{x})
        xl{x}=get(gca,'xlim');
        yl{x}=get(gca,'ylim');
        title(['Neural channel ' ch_str_array{x}])
    end

    linkaxes(ax,'x');
    set(gca,'xlim',[0 time(end)])

    filt_str=input('Apply 4.0 kHz lowpass filter? [y/n], or a vector of channels to filter - [x x x] ','s');
    if strcmp(filt_str,'y')
        low_filt=ones(1,length(ch));
    elseif strcmp(filt_str,'n')
        low_filt=zeros(1,length(ch));
    else
        low_filt=zeros(1,length(ch));
        filt_id=str2num(filt_str);
        for x=1:length(ch)
            if sum(filt_id==ch(x))
                low_filt(x)=1;
            end
        end
    end


    for x=1:length(ch)
        if low_filt(x)
            disp(['Low-pass filtering channel ' num2str(ch(x)) ' to remove movement artifact']);
            neural{x}=bandpass(neural{x},Fs,1,4000, 'hanningfir');
            ax(x+1)=subplot(nfigs,1,x+1);hold on
            plot(time,neural{x},'r')
        end
    end

end

% same command for batchmode or not
if ~batchmode
    if analyze  % compute signal to noise ratio
        if length(ch)>1
            ch_tmp=str2num(input('Which channel(s) to use for data out? For mult channels, syntax is [x x x] ','s'));
            ch_select=intersect(ch_tmp,ch);
            [tmp,id,tmp2]=intersect(ch,ch_tmp);
            invert=invert(id);
            low_filt=low_filt(id);
            for x=1:length(ch)
                z(x)=sum(ch(x)==ch_tmp);
            end
            ch_select=find(z);
            %        ch_select=find(ch==str2num(input('Which channel(s) to use for data out? For mult channels, syntax is [x x x] ','s')))
        else
            ch_select=1;
        end
    else
        ch_select=1;
    end
end

for x=1:length(ch_select)
    if analyze
        if ~batchmode
            in=input(['Threshold for channel ' num2str(ch(ch_select(x)))  '? (hit return for default = ' num2str(TH_default) ' )']);
            if isempty(in)
                TH(x)=TH_default;
            else
                TH(x)=in;
                TH_default=in;
            end
        end
    else
        TH(x)=TH_default;
    end
    id=find(neural{ch_select(x)}>TH(x));
    sdx=diff(sign(diff(neural{ch_select(x)})));
    id2=find(sdx<0)+1;
    %     id_peaks{x}=intersect(id,id2);
    %     peaks{x}=neural{ch_select(x)}(id_peaks{x});
    %     t_peaks{x}=time(id_peaks{x});

    first_der=diff(neural{ch_select(x)});
    not_constant=find(first_der);
    id_peaks{x}=intersect(id,id2);

    % added 12/15/2006 - neural signal must be changing - some problems
    % with signals getting maxxed out and flattened by data acq cutoff.
    %
    % if there are N points in a peak with identical value, this will assign
    % spiketime to the Nth value
    id_peaks{x}=intersect(id_peaks{x},not_constant);
    peaks{x}=neural{ch_select(x)}(id_peaks{x});
    t_peaks{x}=time(id_peaks{x});

    %     if 1
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         % new module added 12/15/2006
    %         needs_reiteration=1;
    %         n_reiterations=0;
    %         % screening out times where >1consec sample is a t_peak
    %         % this will happen, for example, when spikes are truncated by
    %         % data acq magnitude limits.
    %         if 0
    %         find(( t_peaks{x}>12.888&  t_peaks{x}<12.9))
    %                  tmp=t_peaks{x}( t_peaks{x}>12.888&  t_peaks{x}<12.9)
    %                  diff(tmp)*32000
    %
    %                  figure(100);clf;hold on
    %                  neur=neural{ch_select(x)};
    %                  id_tmp=id_peaks{x}(828);
    %                  ax(1)=subplot(3,1,1)
    %                  plot(neur(id_tmp-32:id_tmp+32),'-o')
    %                  ax(2)=subplot(3,1,2)
    %                  plot(sdx(id_tmp-32:id_tmp+32),'r')
    %                  ax(3)=subplot(3,1,3)
    %                  plot(first_der(id_tmp-32:id_tmp+32),'r')
    %                  linkaxes(ax,'x')
    %                  return
    %         end
    %                  while needs_reiteration
    %             % if the difference between two identified peaks is one sample
    %             % (cant do ==1/Fs because of MATLAB precision limit)
    %             temp=find((diff(t_peaks{x})-1/Fs)<1e-6);
    %             if ~isempty(temp);disp(['Found ' num2str(length(temp)) ' double peaks']);end
    % %                  tmp=  [neural{ch_select(x)}(id_peaks{x}(temp)-1) neural{ch_select(x)}(id_peaks{x}(temp)) neural{ch_select(x)}(id_peaks{x}(temp)+1) neural{ch_select(x)}(id_peaks{x}(temp)+2)];
    % %                  sortrows(tmp,1);
    %             if ~isempty(temp)
    %                 for z=1:length(temp)
    % %                    t_peaks{x}(temp(z)+1)=9999999;
    %                     t_peaks{x}(temp(z))=9999999;
    %                     % check to be sure that double peaks are caused by
    %                     % maxed-out (or minned-out) value n
    %                     if abs(neural{ch_select(x)}(id_peaks{x}(temp(z)+1)))~=32768;
    %                         if 0
    %                         figure(100);clf;hold on
    %                         neur=neural{ch_select(x)};
    %                         id_tmp=id_peaks{x}(temp(z));
    %                         plot(neur(id_tmp-32:id_tmp+32),'-o')
    %                         error(['Double peak (' num2str(neural{ch_select(x)}(id_peaks{x}(temp(z)+1))) ') not a maxed-out value - see fig 100']);
    %                         end
    %                         disp(['NOTE: Double peak (' num2str(neural{ch_select(x)}(id_peaks{x}(temp(z)+1))) ') not a maxed-out value']);
    %                     end
    %                 end
    %                 t_peaks{x}=t_peaks{x}(find(t_peaks{x}~=9999999));
    %                  t_peaks{x}( t_peaks{x}>12.888 & t_peaks{x}<12.9)
    %                 n_reiterations=n_reiterations+1;
    %             else
    %                 needs_reiteration=0;
    %             end
    %         end
    %         if n_reiterations
    %             disp(['Performed ' num2str(n_reiterations) ' iterations to get rid of double peaks'])
    %         end
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     end

    if ~batchmode
        subplot(nfigs,1,ch_select(x)+1);hold on
        tbin_for_hist=.05;
        st_quants=floor( t_peaks{x}/tbin_for_hist)+1;
        quantbins=unique(st_quants);
        [b,c]=histc(st_quants,quantbins);
        VEC=sparse(1,ceil(max(time)/tbin_for_hist));
        VEC(quantbins)=b;
        VEC=VEC/max(VEC);   % normalize
        yl_use=yl{ch_select(x)};
        xl_use=xl{ch_select(x)};
        VEC=VEC*.4*diff(yl_use);    % stretch to .4 of the y range
        VEC=VEC+yl_use(1);
        plot((0:(length(VEC)-1))*tbin_for_hist,VEC,'g','linew',2);
        set(gca,'xlim',xl_use);
        set(gca,'ylim',yl_use);
        pause(.2)
    end
end



if analyze
    for z=1:length(ch_select)
        tp= t_peaks{z};
        ct=1;
        id=[];
        tic
        id_ones_near_spike=zeros(size(time));
        for x=round(-.001*Fs):round(.001*Fs)
            id_tmp=id_peaks{z}+x;
            id_tmp=nonzeros((id_tmp+abs(id_tmp))/2); % clear times before beginning
            id_tmp=id_tmp(find(id_tmp<=length(time)));       % clear times after end
            id_ones_near_spike(id_tmp)=1;
        end
        id_ones_near_spike=find(id_ones_near_spike);
        id=id_ones_near_spike;
        neural_spikes=peaks{z};
        neural_nospikes=neural{ch_select(z)}(setdiff(1:length(time),id));
        mean_spikes=mean(neural_spikes);
        mean_nospikes=mean(neural_nospikes);
        std_spikes=std(neural_spikes);
        std_nospikes=std(neural_nospikes);
        spikes_mag=mean_spikes;
        nospikes_mag=mean_nospikes+2*std_nospikes;
        SNR(z)=spikes_mag/nospikes_mag;
        disp(['Signal to noise ratio for channel ' num2str(ch(ch_select(z))) ' = ' num2str(SNR(z))])
        if ~batchmode
            subplot(nfigs,1,ch_select(z)+1);hold on
            xl=get(gca,'xlim');
            yl=get(gca,'ylim');
            disp(['Signal to noise ratio for channel ' num2str(ch(ch_select(z))) ' = ' num2str(SNR(z))])
            plot(xl,TH(z)*ones(1,2),'k--')
            plot(xl,spikes_mag*ones(1,2),'c')
            plot(xl,nospikes_mag*ones(1,2),'y')
            plot(xl,mean_nospikes*ones(1,2),'y:')
        end
    end
end

if analyze
    Data.file=file;
    Data.TH=TH;
    Data.Fs=Fs;
    Data.total_n_samples=length(neural{1});
    Data.channels=ch(ch_select);
    Data.SNR=SNR;
    Data.inverted_waveform=invert;
    Data.lowpass_filter=low_filt;
    for x=1:length(ch_select)
        Data.spiketimes{ch(ch_select(x))}=t_peaks{x};
        if ~batchmode
            in=input(['Is channel ' num2str(ch(ch_select(x))) ' a single unit?  (y/n)'],'s');
            if in=='y'
                Data.is_single_unit(x)=1;
            elseif in=='n'
                Data.is_single_unit(x)=0;
            else
                error
            end
            Data.analyzed_in_batchmode=0;
        else
            Data.is_single_unit(x)=9999;
            Data.analyzed_in_batchmode=1;
        end
    end
    eval(sprintf('save %s.neuralnot.mat Data',file))
end
