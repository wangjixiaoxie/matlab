load UDsaved021011y50g89.mat
%%%%%%
winwidth=400;
windowstarts=-4000:100:4000;
for thissyllable=1:3
    for mj=1:length(windowstarts)
        windowstart=windowstarts(mj)
        windowend=windowstart+winwidth;
        clear alldata
        numrends=zeros(1,12);
        for thisfolder=1:12
            count=zeros(1,40);
            % What is the range of spikes?
            for goodclust=clusters{thisfolder*2}
                clear pdata
                thisclust=goodclust;
                for thissong=1:length(neuron(thisfolder).song)
                    % if this song has all clusters:
                    if thisclust<length(neuron(thisfolder).song(thissong).spiketimes)+1
                        if length(neuron(thisfolder).song(thissong).pitchdata)>thissyllable-1
                            for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                numrends(thisfolder)=numrends(thisfolder)+1;
                                % pitch contour
                                pdata=(neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:));
                                thesespikes=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                    -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                spikes=length(thesespikes(find(thesespikes>windowstart & thesespikes<windowend)));
                                count(spikes+1)=count(spikes+1)+1;
                                alldata{thisfolder,spikes+1}.datt(count(spikes+1),:)=pdata;
                            end
                        end
                    end
                end
            end
        end
        % What is distn of spike rates for each neuron
        clear thesesizes
        clear foldermean
        for thisfolder=1:12
            for j=1:size(alldata,2)
                if isempty(alldata{thisfolder,j})
                    thesesizes(j)=0;
                else
                    thesesizes(j)=size(alldata{thisfolder,j}.datt,1);
                end
            end
            % foldermean(thisfolder)=dot([0:1:size(alldata,2)-1],thesesizes)/sum(thesesizes);
            foldermean(thisfolder)=dot([1:1:size(alldata,2)-1],thesesizes(2:end))/sum(thesesizes(2:end));            
        end
        %
        pitchesDiff=[];
        for thisfolder=1:12
            pitchdataLow=[];
            pitchdataHigh=[];
            for i=2:size(alldata,2) % ignore zero's
%             for i=1:size(alldata,2)
                if ~isempty(alldata{thisfolder,i})
                    if i-1<foldermean(thisfolder)
                        pitchdataLow=[pitchdataLow;alldata{thisfolder,i}.datt];
                    else
                        pitchdataHigh=[pitchdataHigh;alldata{thisfolder,i}.datt];
                    end
                end
            end
            pitchesDiff(thisfolder,:)=mean(pitchdataHigh)-mean(pitchdataLow);
        end
        weightedavg(mj,thissyllable).datt=(numrends*pitchesDiff)/sum(numrends); % weights each neuron by the number of renditions (seems reasonable)
    end
end

figure;hold on;
mj=1
for i=1:size(weightedavg,1)
    mn1s(i)=mean(weightedavg(i,mj).datt(250:400));
end
plot(windowstarts+winwidth/2-250/8,mn1s,'b')
mj=2
for i=1:size(weightedavg,1)
    mn2s(i)=mean(weightedavg(i,mj).datt(600:1000));
end
plot(windowstarts+winwidth/2-600/8,mn2s,'r')
mj=3
for i=1:size(weightedavg,1)
    mn3s(i)=mean(weightedavg(i,mj).datt(700:850));
end
plot(windowstarts+winwidth/2-700/8,mn3s,'k')
mj=3
for i=1:size(weightedavg,1)
    mn4s(i)=mean(weightedavg(i,mj).datt(1350:1450));
end
plot(windowstarts+winwidth/2-1350/8,mn4s,'g')







% Song correlations analysis
clear all
load UDsaved021011y50g89.mat
neuron(12).song=[neuron(12).song(1:11) neuron(12).song(13:14)]; % removes song with no labels
% Part 1: Pitch correlations
numsyllables=3;
        % A. Calculate mean pitch of each syllable irregardless of where syllable occurs in song and store in mnpd  
                % Where to look in each syllable for the harmonic stack
                        clear ptwins
                        ptwins(1).datt=250:400;
                        ptwins(2).datt=600:1000;
                        ptwins(3).datt=700:850;
                % Creates variable pdata, containing mean pitch of each rendition
                % (columns) of each syllable (rows)
                        clear pdata
                        for thissyllable=1:numsyllables % each syllable
                            count=0;
                            for thisfolder=1:length(neuron) % each neuron
                                for thissong=1:length(neuron(thisfolder).song) % each song
                                    if length(neuron(thisfolder).song(thissong).pitchdata)>thissyllable-1 % if the song has each syllable
                                        for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable}) % for each syll rendition
                                            count=count+1; % increment
                                            pdata(thissyllable,count)=mean(neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt...
                                                (thisrendition,ptwins(thissyllable).datt)); % mean pitch of harmonic stack
                                        end
                                    end
                                end
                            end
                        end
                % Creates variable mnpd, containing mean pitch for each syllable
                    for i=1:numsyllables;mnpd(i)=mean(pdata(i,1:length(find(pdata(i,:)>0))));end
        % B. Create FFcorrs
        %   The structure FFcorrs contains a "data" field for each song,
        % containing time relative to first syllable of that song,FF,and 
        % syllable identity for each rendition.
        %   FFcorrs allows calculations of FF vs. time in song, respective or
        % irrespective of syllable identity.
                    FFcorrs=[];
                    songcount=0;
                    for thisfolder=1:length(neuron) % for each folder
                        for thissong=1:length(neuron(thisfolder).song) % for each song
                            songcount=songcount+1;
                            numcount=0;
                            %%%%%%% when is the time of the earliest syllable in this song (all times for this song will be relative to this)
                                    allons=[];
                                    for kk=1:length(neuron(thisfolder).song(thissong).noteons)
                                        allons=[allons neuron(thisfolder).song(thissong).noteons{1}()];
                                    end
                                    mntimes=min(allons); % amount of time from song onet to first syllable onset (to be subtracted from time variables for all syllables)
                            %%%%%%%
                            for thissyllable=1:length(neuron(thisfolder).song(thissong).noteons) % for each syllable
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable}) % for each rendition
                                    numcount=numcount+1; % each rendition
                                    tshift=mean(ptwins(thissyllable).datt/8)+16; % amount of time (ms) from syllable onset to harmonic stack (to be added onto time variable)
                                    % 1=times of harmonic stack onset relative to onset of first syllable in this song
                                    FFcorrs(songcount).data(1,numcount)=tshift-mntimes+neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition); % times
                                    % 2=FF for that rendition (as difference from mean for all renditions of that syllable)
                                    FFcorrs(songcount).data(2,numcount)=mean(neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,ptwins(thissyllable).datt))...
                                        -mnpd(thissyllable);
                                    % 3=Syllable number (e.g. 1 corresponds to A, 2 corresponds to B, etc.)
                                    FFcorrs(songcount).data(3,numcount)=thissyllable;
                                end
                            end
                        end
                    end

        % C. Plot time in song vs. FF deviation ---- THERE IS A DOWNWARD TREND
                    figure;hold on;count=0;for i=1:length(FFcorrs);[a,b]=sort(FFcorrs(i).data(1,:));plot(FFcorrs(i).data(1,b),FFcorrs(i).data(2,b),'.');end
% Part 2: SPIKE correlations
        % A. Calculate mean number of spikes around each syllable, irrespective of where that syllable occurs in song.
        % I do not pool across neurons (i.e. folders), thus no weightedaverage
                winstart=-200; % Start of window (relative to harmonic stack onset) where I look for spikes for that neuron
                winend=0; % End of window (relative to harmonic stack onset) where I look for spikes for that neuron
                for thissyllable=1:numsyllables % for each syllable
                    clear foldermn folderN
                    for thisfolder=1:length(neuron) % for each neuron
                        goodclust=clusters{thisfolder*2};
                        spikecounts=[];
                        count=0;
                        % What is the range of spikes?
                        for thissong=1:length(neuron(thisfolder).song) % for each song
                            if length(neuron(thisfolder).song(thissong).pitchdata)>thissyllable-1 % if there are renditions of the syllable in this song
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable}) % for each rendition
                                    tshift=mean(ptwins(thissyllable).datt/8)+16; % amount of time (ms) from syllable onset to harmonic stack (to be added onto time variable)
                                    % Count the spikes in the window
                                    spikes=0;
                                    for jk=1:length(goodclust)
                                        thisclust=goodclust(jk);
                                        % if this song has all clusters:
                                        if thisclust<length(neuron(thisfolder).song(thissong).spiketimes)+1
                                            thesespikes=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                                -(tshift+neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition)); % all spikes in the song - times relative to syllable onset
                                            spikes=spikes+length(thesespikes(find(thesespikes>winstart & thesespikes<winend)));
                                        end
                                    end
                                    count=count+1;
                                    spikecounts(count)=spikes;
                                end
                            end
                        end
                        foldermn(thisfolder)=mean(spikecounts);
                        folderN(thisfolder)=length(spikecounts);
                    end
                    foldermean(thissyllable,:)=foldermn; % average
                end
            % B. mean spike count for each syllable and each neuron/folder
                    SPKcorrs=[];
                    songcount=0;
                    for thisfolder=1:length(neuron) % for each folder
                        goodclust=clusters{thisfolder*2};
                        for thissong=1:length(neuron(thisfolder).song) % for each song
                            numcount=0;
                            songcount=songcount+1;
                            %%%%%%% when is the time of the earliest syllable in this song (all times for this song will be relative to this)
                                    allons=[];
                                    for kk=1:length(neuron(thisfolder).song(thissong).noteons)
                                        allons=[allons neuron(thisfolder).song(thissong).noteons{1}()];
                                    end
                                    mntimes=min(allons);
                            %%%%
                            for thissyllable=1:length(neuron(thisfolder).song(thissong).noteons) % for each syllable
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable}) % for each rendition
                                    tshift=mean(ptwins(thissyllable).datt/8)+16; % amount of time (ms) from syllable onset to harmonic stack (to be added onto time variable) 
                                    % 1=times of harmonic stack onset relative to onset of first syllable in this song
                                    numcount=numcount+1;
                                    SPKcorrs(songcount).data(1,numcount)=tshift+neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition)-mntimes;
                                    % Count the spikes in the window
                                    spikes=0;
                                    for jk=1:length(goodclust)
                                        thisclust=goodclust(jk);
                                        % if this song has all clusters:
                                        if thisclust<length(neuron(thisfolder).song(thissong).spiketimes)+1
                                            thesespikes=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                                -(tshift+neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition)); % all spikes in the song - times relative to syllable onset
                                            spikes=spikes+length(thesespikes(find(thesespikes>winstart & thesespikes<winend)));
                                        end
                                    end
                                    % 2=Spikes at that time (as deviation from mean for that neuron and that syllable)
                                    SPKcorrs(songcount).data(2,numcount)=spikes-foldermean(thissyllable,thisfolder);%%%%mean(neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,ptwins(thissyllable).datt))-mnpd(thissyllable);
                                    % 3=Syllable identity (e.g. 1=syll A)
                                    SPKcorrs(songcount).data(3,numcount)=thissyllable;
                                end
                            end
                        end
                    end

figure;hold on;
%%%%%%%%%%% Plot all data
    subplot(311);hold on;
    count=0;
    alltimes=[];
    allFFs=[];
    allspks=[];
    for i=1:length(SPKcorrs)
        [a,b]=sort(SPKcorrs(i).data(1,:));
        plot(SPKcorrs(i).data(1,b),SPKcorrs(i).data(2,b),'.')
        alltimes=[alltimes FFcorrs(i).data(1,:)];
        allFFs=[allFFs FFcorrs(i).data(2,:)];
        allspks=[allspks SPKcorrs(i).data(2,:)];
    end
    subplot(312);hold on;count=0;for i=1:length(FFcorrs);[a,b]=sort(FFcorrs(i).data(1,:));plot(FFcorrs(i).data(1,b),FFcorrs(i).data(2,b),'.');end
%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(313);hold on;
winwidthms=300;
winstepms=50;
clear mnspk mntm mnFF
for i=1:floor(2000/winstepms)
    mnspk(i)=mean(allspks(find(alltimes>(i-1)*winstepms & alltimes < (i-1)*winstepms+winwidthms)));
    mntm(i)=mean(alltimes(find(alltimes>(i-1)*winstepms & alltimes < (i-1)*winstepms+winwidthms)));
    mnFF(i)=mean(allFFs(find(alltimes>(i-1)*winstepms & alltimes < (i-1)*winstepms+winwidthms)));
    mnN(i)=length(find(alltimes>(i-1)*winstepms & alltimes < (i-1)*winstepms+winwidthms));
end
plot(mntm,mnFF,'b');hold on;plot(mntm,50*mnspk,'r')


figure;hold on;plot(mntm,mnFF,'b');hold on;plot(mntm,50*mnspk,'r')
figure;plot(xcov(mnspk,mnFF))
% [a,b]=sort(alltimes);
% figure;hold on;plot(runningaverage(alltimes(b),50),50*(runningaverage(allFFs(b),50)-1),'b')
% plot(runningaverage(alltimes(b),50),runningaverage(allspks(b),50),'r')
% figure;plot(xcov(runningaverage(allFFs(b(1:200)),50),runningmedian(allspks(b(1:200)),50)))
