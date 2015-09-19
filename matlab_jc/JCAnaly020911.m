
% .song.mat
    % Song data, sampled as 32kHz
% .spk
    % First column = cluster number (for each neuron, only a subset of
            % clusters are actual spikes - and the rest are just noise -
            % see e-mail from raghav)
    % Second column = spike time (seconds)
    
%     y50g89_081610160745.song.mat





% Undirected
load /cardinal/HVCAnalysedData/AllSavedFeb9.mat

% Directed


basedir='/cardinal/HVCAnalysedData/'
a=1;folders{a}=[basedir 'y50g89_081610_4/Directed'];clusters{a}=[1 3 4 5 6 7];prefixes{a}='Chan3_';
a=2;folders{a}=[basedir 'y50g89_081610_4/UnDirected'];clusters{a}=[1 3 4 5 6 7];prefixes{a}='Chan3_';
a=3;folders{a}=[basedir 'y50g89_081610_5/Directed'];clusters{a}=[1 2];prefixes{a}='Chan3_';
a=4;folders{a}=[basedir 'y50g89_081610_5/UnDirected'];clusters{a}=[1 2];prefixes{a}='Chan3_';
a=5;folders{a}=[basedir 'y50g89_081710_6/Directed'];clusters{a}=[3 4 5];prefixes{a}='Chan2_';
a=6;folders{a}=[basedir 'y50g89_081710_6/UnDirected'];clusters{a}=[3 4 5];prefixes{a}='Chan2_';
a=7;folders{a}=[basedir 'y50g89_081710_7/Directed'];clusters{a}=[2];prefixes{a}='Chan2_';
a=8;folders{a}=[basedir 'y50g89_081710_7/UnDirected'];clusters{a}=[2];prefixes{a}='Chan2_';
a=9;folders{a}=[basedir 'y50g89_081710_8/Directed'];clusters{a}=[1 6 9 12];prefixes{a}='Chan2_';
a=10;folders{a}=[basedir 'y50g89_081710_8/UnDirected'];clusters{a}=[1 6 9 12];prefixes{a}='Chan2_';
a=11;folders{a}=[basedir 'y50g89_081810_9/Directed'];clusters{a}=[3];prefixes{a}='Chan3_';
a=12;folders{a}=[basedir 'y50g89_081810_9/UnDirected'];clusters{a}=[3];prefixes{a}='Chan3_';
a=13;folders{a}=[basedir 'y50g89_082010_10/Directed'];clusters{a}=[4 7 9];prefixes{a}='Chan3_';
a=14;folders{a}=[basedir 'y50g89_082010_10/UnDirected'];clusters{a}=[4 7 9];prefixes{a}='Chan3_';
a=15;folders{a}=[basedir 'y50g89_082010_11/Directed'];clusters{a}=[1 6 9 10];prefixes{a}='Chan3_';
a=16;folders{a}=[basedir 'y50g89_082010_11/UnDirected'];clusters{a}=[1 6 9 10];prefixes{a}='Chan3_';
a=17;folders{a}=[basedir 'y50g89_082010_12/Directed'];clusters{a}=[2 4 7 9 12];prefixes{a}='Chan4_';
a=18;folders{a}=[basedir 'y50g89_082010_12/UnDirected'];clusters{a}=[2 4 7 9 12];prefixes{a}='Chan4_';
a=19;folders{a}=[basedir 'y50g89_082310_13/Directed'];clusters{a}=[2 5];prefixes{a}='Chan4_';
a=20;folders{a}=[basedir 'y50g89_082310_13/UnDirected'];clusters{a}=[2 5];prefixes{a}='Chan4_';
a=21;folders{a}=[basedir 'y50g89_082410_14/Directed'];clusters{a}=[5 7 10];prefixes{a}='Chan5_';
a=22;folders{a}=[basedir 'y50g89_082410_14/UnDirected'];clusters{a}=[5 7 10];prefixes{a}='Chan5_';
a=23;folders{a}=[basedir 'y50g89_101210_15/Directed'];clusters{a}=[1 3];prefixes{a}='Chan4_';
a=24;folders{a}=[basedir 'y50g89_101210_15/UnDirected'];clusters{a}=[1 3];prefixes{a}='Chan4_';



for mm=1:2:23%:10%each folder
    clearvars -except  mm folders neuron clusters prefixes allvals
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
    n=(m+1)/2;
    neuron(n).song=song;
end

% neuron(folder)
%	 .song(song)
        % .Song --> entire song file
        % .noteons{syllable}(rendition)   % note onsets (minus 30ms) relative to song onset (30ms subtracted to allow pitch quant)
        % .noteoffs{syllable}(rendition)  % note offsets (minus 30ms) relative to song onset
        % .spiketimes{cluster}            % spike times (minus 30ms) relative to song onset (30ms subtracted to align with note onsets)
        % .notedata{syllable}.datt(rendition,:)   % raw data for note - used to calculate pitch contours
        % .pitchdata{syllable}.datt(rendition,:)  % pitch contours - first point is 16ms after note onset (512 points)
%




% PLOT RRRRRasters!
pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
vals=1300;
figure;hold on
for thisfolder=1:11%length(neuron)
    goodclust=clusters{thisfolder*2};
    for thissyllable=1
        for thissong=1:length(neuron(thisfolder).song)
            plot(pitchtimes,neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt)
            for thisclust=1:length(goodclust)
                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                thesespikes=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                    -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                if ~isempty(thesespikes)
                    vals=vals-1;
                plot(thesespikes,vals,'b.')
                end
                end
            end
        end
    end
end



wins=[-820:10:820];

for mls=1:length(wins);
    winstart=wins(mls);
    winstop=wins(mls)+50;

%%%%% SYLLABLE A
                % Analyze - this loop returns:
                    % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
                    % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
                % for the first syllable, the harmonic part goes from 45-70ms (includes 16ms adjustment - so really 29-54ms)
                pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
                harmstart=45;
                harmend=70;
                neuroncount=0;
                clear spikecount pitchmean 
                for thisfolder=1:11 % for the neuron/recording wire
                    goodclust=clusters{thisfolder*2};
                    for thissyllable=1 % for the syllable
                        for thisclust=1:length(goodclust)     % for the cluster
                            neuroncount=neuroncount+1; % each unit
                            count=0; % each rendition with this unit
                            for thissong=1:length(neuron(thisfolder).song)
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                    timewindow=harmstart+winstart:harmend+winstop;
                                    spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                        -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                    spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end

                % Only analyze results for the neurons that actually spike in the vicinity
                % of the harmonic stack note - the variable 'doesspike' evaluates this
                clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=sum(spikecount{neuralunit}())>10; % nine units
                    end

                    % Calculate correlation coefficients and corresponding t-stats
                    clear rval tval
                        for neuralunit=1:length(spikecount) % for each 'neuron'
                            if doesspike(neuralunit)
                                m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                                r=m(2);
                                rval(neuralunit)=r;
                                tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                            end
                        end
                        rval(doesspike)
                        tval(doesspike)
                SyllA.rval=rval;
                SyllA.tval=tval;
                SyllA.spikecount=spikecount;
                SyllA.pitchmean=pitchmean;
                
%%%%% SYLLABLE B
                % Analyze - this loop returns:
                    % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
                    % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
                % for the second syllable, the harmonic part goes from 110-140ms (includes 16ms adjustment - so really 94-124ms)
                clear spikecount pitchmean
                pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
                harmstart=80; % quite a bit of leeway here
                harmend=160; % quite a bit of leeway here
                neuroncount=0;
                for thisfolder=1:11 % for the neuron/recording wire
                    goodclust=clusters{thisfolder*2};
                    for thissyllable=1 % for the syllable
                        for thisclust=1:length(goodclust)     % for the cluster
                            neuroncount=neuroncount+1; % each unit
                            count=0; % each rendition with this unit
                            for thissong=1:length(neuron(thisfolder).song)
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                    timewindow=harmstart+winstart:harmend+winstop;
                                    spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                        -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                    spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end

                % Only analyze results for the neurons that actually spike in the vicinity
                % of the harmonic stack note - the variable 'doesspike' evaluates this
                clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=sum(spikecount{neuralunit}())>10; % nine units
                    end
                    clear rval tval
                    % Calculate correlation coefficients and corresponding t-stats
                        for neuralunit=1:length(spikecount) % for each 'neuron'
                            if doesspike(neuralunit)
                                m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                                r=m(2);
                                rval(neuralunit)=r;
                                tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                            end
                        end
                        rval(doesspike)
                        tval(doesspike)
                SyllB.rval=rval;
                SyllB.tval=tval;
                SyllB.spikecount=spikecount;
                SyllB.pitchmean=pitchmean;
    
    
 %%%%% SYLLABLE C
                % Analyze - this loop returns:
                    % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
                    % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
                % for the first syllable, the harmonic part goes from 45-70ms (includes 16ms adjustment - so really 29-54ms)
                pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
                harmstart=90;
                harmend=130;
                neuroncount=0;
                clear spikecount pitchmean     
                for thisfolder=1:11 % for the neuron/recording wire
                    goodclust=clusters{thisfolder*2};
                    for thissyllable=1 % for the syllable
                        for thisclust=1:length(goodclust)     % for the cluster
                            neuroncount=neuroncount+1; % each unit
                            count=0; % each rendition with this unit
                            for thissong=1:length(neuron(thisfolder).song)
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                    timewindow=harmstart+winstart:harmend+winstop;
                                    spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                        -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                    spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end

                % Only analyze results for the neurons that actually spike in the vicinity
                % of the harmonic stack note - the variable 'doesspike' evaluates this
                clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=sum(spikecount{neuralunit}())>10; % nine units
                    end

                    % Calculate correlation coefficients and corresponding t-stats
                    clear rval tval
                        for neuralunit=1:length(spikecount) % for each 'neuron'
                            if doesspike(neuralunit)
                                m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                                r=m(2);
                                rval(neuralunit)=r;
                                tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                            end
                        end
                        rval(doesspike)
                        tval(doesspike)
                SyllC.rval=rval;
                SyllC.tval=tval;
                SyllC.spikecount=spikecount;
                SyllC.pitchmean=pitchmean;
   
mnabsA(mls)=median(abs(SyllA.tval(SyllA.tval~=0 & ~isnan(SyllA.tval))));
mnabsB(mls)=median(abs(SyllB.tval(SyllB.tval~=0 & ~isnan(SyllB.tval))));
mnabsC(mls)=median(abs(SyllC.tval(SyllC.tval~=0 & ~isnan(SyllC.tval))));
mnA(mls)=median((SyllA.tval(SyllA.tval~=0 & ~isnan(SyllA.tval))));
mnB(mls)=median((SyllB.tval(SyllB.tval~=0 & ~isnan(SyllB.tval))));
mnC(mls)=median((SyllC.tval(SyllC.tval~=0 & ~isnan(SyllC.tval))));
seA(mls)=std((SyllA.tval(SyllA.tval~=0 & ~isnan(SyllA.tval))))/sqrt(sum(SyllA.tval~=0 & ~isnan(SyllA.tval)));
seB(mls)=std((SyllB.tval(SyllB.tval~=0 & ~isnan(SyllB.tval))))/sqrt(sum(SyllB.tval~=0 & ~isnan(SyllB.tval)));
seC(mls)=std((SyllC.tval(SyllC.tval~=0 & ~isnan(SyllC.tval))))/sqrt(sum(SyllC.tval~=0 & ~isnan(SyllC.tval)));


end
    

figure;plot([-1000:10:1000],zeros(1,201),'k-')
hold on;plot(wins,mnA-seA)
hold on;plot(wins,mnA+seA)
hold on;plot(wins,mnB+seB,'g')
hold on;plot(wins,mnB-seB,'g')
hold on;plot(wins,mnC+seC,'r')
hold on;plot(wins,mnC-seC,'r')