function [neuron]=HVCprocessing_y50g89b(folders,clusters,prefixes)
% neuron(folder)
%	 .song(song)
        % .Song --> entire song file
        % .noteons{syllable}(rendition)   % note onsets (minus 30ms) relative to song onset (30ms subtracted to allow pitch quant)
        % .noteoffs{syllable}(rendition)  % note offsets (minus 30ms) relative to song onset
        % .spiketimes{cluster}            % spike times (minus 30ms) relative to song onset (30ms subtracted to align with note onsets)
        % .notedata{syllable}.datt(rendition,:)   % raw data for note - used to calculate pitch contours
        % .pitchdata{syllable}.datt(rendition,:)  % pitch contours - first point is 16ms after note onset (512 points)


mmat=[1:1:24];        
for mm=mmat %each folder
    clearvars -except  mm folders neuron clusters prefixes allvals isDir
    m=mm;
    currfolder=folders{m};
    cmd=['cd ' currfolder]
    eval(cmd);
    goodclusts=clusters{m};
    chanprefix=prefixes{m};
    prewinsong=30; % ms,for pitch quant
    prewinspk=50; % add on prewinsong - so really 80ms
    postwinspk=50;
    ntlength=8000;
    dirf('*.song.mat','batchsongs')
    NFFT=1020;
    note_cnt=0;avnote_cnt=0;fcnt=0;
    ff=load_batchf('batchsongs');
    eachind=0;
    for ifn=1:length(ff) % for each song
        fn=ff(ifn).name(1:end-9);
        fnn=[fn,'.not.mat'];
        if (~exist(fnn,'file'))
            continue;
        end
        load(fnn)
        labels = lower(labels);
        labels(findstr(labels,'0'))='-';
        cmd=['load ' ff(ifn).name];
        eval(cmd);
        % spikes
        fid=fopen([chanprefix fn '.spk']); % change chan prefix for each bird
        flsv=[];
        count=0;
        while (1)
            if (length(flsv)==0)
                fl = fgetl(fid);
            else
                fl = flsv;
                flsv=[];
            end
            if (~ischar(fl))
                break;
            end
            if length(fl)==0
                continue;
            end
            count=count+1;
            placeholder=find(isspace(fl));
            ind_clust(count)=str2num(fl(1:placeholder-1)); % cluster number
            spiketm(count)=str2num(fl(placeholder+1:end)); % spike time (seconds)
        end
        
        % This loop pulls out the spike times for each cluster (in the
        % song that is currently open) --> prewinsong is subtracted
        for jj=1:max(ind_clust)
            song(ifn).spiketimes{jj}=spiketm(ind_clust==jj)*1e3-prewinsong;
        end
        
        % This takes the raw song
        dat=Song;
        song(ifn).Song=dat;
        
        % This gets the onsets/offsets for the structure and processes data
        % for later pitch measurement
        allnts='abc';%unique(labels);
        windowingFF(1,:)=[1300 1500];windowingFF(2,:)=[1350 1550];windowingFF(3,:)=[1300 1500];     
        for i=1:length(allnts) % for each note present in this song
            currnt=allnts(i);
            % excludes notes that start too early or late in song (first 80ms)
            ind1=find(onsets>prewinsong & offsets*1e-3+0.25<length(Song)/Fs);
            onsets1=onsets(ind1);
            offsets1=offsets(ind1);
            p=findstr(labels(ind1),currnt);
            for ii=1:length(p) % for each instance of that note in this song
                ton=onsets1(p(ii))-prewinsong; % starting 30ms before syllable onset
                toff=offsets1(p(ii));
                ti1=ceil((ton*1e-3)*Fs);
                ti2=ti1+ntlength;
                datt=dat(ti1:ti2);
                song(ifn).noteons{i}(ii)=ton;
                song(ifn).noteoffs{i}(ii)=toff;  
                song(ifn).notedata{i}.datt(ii,:)=datt;                  
              % calculate spectrogram 
                    N=1024;SAMPLING=Fs;sigma=64;OVERLAP=1020;
                    %sonodvA program
                    t=-N/2+1:N/2;
                    sigma=(sigma/1000)*SAMPLING;
                    %Gaussian and first derivative as windows.
                    w=exp(-(t/sigma).^2);
                    timebins=floor((length(datt)-(N))/(N-OVERLAP))+1;
                    freqbins=(N/2)+1;
                    Nyquist=SAMPLING/2;
                    step=Nyquist/freqbins; %This is the width (in ms) between bins.
                    songram=zeros(freqbins,timebins);
                    songram=abs(flipdim(spectrogram(datt,w,OVERLAP,N,SAMPLING),1)); %gaussian windowed spectrogram
              % calculate pitch     
                    mini=round(windowingFF(i,2)/step);
                    maxi=round(windowingFF(i,1)/step);
                    % for each rendition of that syllable
                    harms=[1 2 3];%floor(Nyquist/windowingFF(ii,2));
                    for currenttime_bin=1:size(songram,2)
                        slice=songram(:,currenttime_bin); %power at each frequency for the current time bin
                        Powerful=[];
                        Freqbinest=[];
                        Powerful(1)=0;
                        Freqbinest(1)=0;
                        for iii=1:length(harms)
                            minimum=freqbins-mini*harms(iii);
                            maximum=freqbins-maxi*harms(iii);
                            freq_window=slice(minimum:maximum);
                            [Pow,Ind]=max(freq_window);
                            % Interpolation--subtraction of 0.5 because we care about the
                            % central frequency of the bin, not the upper boundary.
                            if Ind==1 || Ind==length(freq_window)
                                Indest=Ind;
                            else
                                Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
                            end
                            Real_Index=minimum+Indest-1; % -1 to account for window;
                            Freqbinest(iii)=(freqbins-Real_Index)/harms(iii);
                            Powerful(iii)=Pow;
                        end
                        normalizer=sum(Powerful);
                        normpower=Powerful/normalizer;
                        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
                    end
              song(ifn).pitchdata{i}.datt(ii,:)=freqbin_estimate*(Nyquist/(freqbins-1));    
            end
        end
    end
    if isDir
        n=(m+1)/2;
    else
        n=m/2;
    end
    neuron(n).song=song;
end

