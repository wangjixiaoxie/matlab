% steps for saucy
% label files

% Which notes do I want to analyze?
n=1;
avls(n).birdname='r30g38';
nts='ba'
avls(n).nts=nts;
chan='1' % 0 for most, 1 for bl82bl81
channame='obs1';
SAUCY_beta('r30g38_140411_140412.414.cbin',chan,2) % all default except set the high/low amplitude limits just outside noise
SAUCY_beta('r30g38_140411_140412.414.cbin',chan,'batchmode')
% 


%%% NOTE A
                chan='1'
                fid=fopen('batchNeural');
                nt='a';
                count=0;
                countA=0;
                 [b,a]=butter(4,[300/32000],'high');
                 clear SongA
                 clear spikesArel waveformArel
                 fs=32000;
                    while 1
                        count=count+1;
                        fn=fgetl(fid);
                        if (~ischar(fn));break;end
                        disp(['Loading ' fn '.not.mat and .neuralnot.mat'])
                        eval(sprintf('load %s.not.mat',fn))
                        eval(sprintf('load %s.neuralnot_CH1.mat',fn))
                        SongA(count).onsets=onsets(find(labels==nt));
                        SongA(count).spiketimes=Data.spiketimes{2}*1000;
                        filteredsig=filtfilt(b,a,evsoundin('',fn,channame));
                        SongA(count).rawfiltwaveform=abs(filteredsig-mean(filteredsig));
                        for i=1:length(SongA(count).onsets)
                            countA=countA+1;
                            spikesArel{countA}=SongA(count).spiketimes-SongA(count).onsets(i);
                            waveformArel(countA,:)=SongA(count).rawfiltwaveform(SongA(count).onsets(i)*(fs/1000)-fs*2:SongA(count).onsets(i)*(fs/1000)+fs);
                        end
                    end
                figure;hold on;for i=1:length(spikesArel);plot(spikesArel{i},-50*i,'b.');end;xlim([-2000 1000])


        % is there song-locked activity?
        % RAW FILTERED SIGNAL
                    ravg1=runningaverage(mean(waveformArel(39:76,:)),3200);
                    ravg2=runningaverage(mean(waveformArel(1:38,:)),3200);
                    figure;hold on;plot(ravg1,'b');plot(ravg2,'r')
                    corrcoef(ravg1,ravg2)

        % SIGNAL AFTER SAUCY
                clear timevals firingrateA
                for i=1:length(spikesArel)
                    timewinstart=[-1025:50:975];
                    timewinwidth=100;
                    for j=1:length(timewinstart)
                        firingrateA(i,j)=length(find(spikesArel{i}>timewinstart(j) & spikesArel{i}<timewinstart(j)+timewinwidth));
                        timevals(j)=timewinstart(j)+timewinwidth/2;
                    end
                end
                figure;hold on;
                plot(timevals,median(firingrateA(1:round(length(spikesArel)/2),:)))
                plot(timevals,median(firingrateA(round(length(spikesArel)/2)+1:end,:)),'r')
                corrcoef(median(firingrateA(1:round(length(spikesArel)/2),:)),median(firingrateA(round(length(spikesArel)/2)+1:end,:)))
%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE B %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
                chan='1'
                fid=fopen('batchNeural');
                nt='b';
                count=0;
                countB=0;
                 [b,a]=butter(4,[300/32000],'high');
                 clear SongB 
                 clear spikesBrel waveformBrel
                 fs=32000;
                    while 1
                        count=count+1;
                        fn=fgetl(fid);
                        if (~ischar(fn));break;end
                        disp(['Loading ' fn '.not.mat and .neuralnot.mat'])
                        eval(sprintf('load %s.not.mat',fn))
                        eval(sprintf('load %s.neuralnot_CH1.mat',fn))
                        SongB(count).onsets=onsets(find(labels==nt));
                        SongB(count).spiketimes=Data.spiketimes{2}*1000;
                        filteredsig=filtfilt(b,a,evsoundin('',fn,channame));
                        SongB(count).rawfiltwaveform=abs(filteredsig-mean(filteredsig));
                        for i=1:length(SongB(count).onsets)
                            countB=countB+1;
                            spikesBrel{countB}=SongB(count).spiketimes-SongB(count).onsets(i);
                            waveformBrel(countB,:)=SongB(count).rawfiltwaveform(SongB(count).onsets(i)*(fs/1000)-fs*2:SongB(count).onsets(i)*(fs/1000)+fs);
                        end
                    end
                figure;hold on;
                for i=1:length(spikesBrel);
                    if sum(isHIT2==i)>0
                        plot(spikesBrel{i},-50*i,'r.');xlim([-2000 1000])
                    else
                        plot(spikesBrel{i},-50*i,'b.');xlim([-2000 1000])
                    end
                end
                
       % HIT vs ESCAPE
               pretimems=0;
               fvalsSongB=findwnoteSPK('batchNeuralnotes','b','','',0,[2000 2700],1e4,1,'obs0',0,pretimems);
               clear isTRIG isCATCH
               for i=1:length(fvalsSongB);isTRIG(i)=fvalsSongB(i).TRIG;end
               for i=1:length(fvalsSongB);isCATCH(i)=fvalsSongB(i).CATCH;end
               isHIT2=find(isTRIG & ~isCATCH); % for identifying fvals
               isESC2=find(~isTRIG | isCATCH);

        % is there song-locked activity?
        % RAW FILTERED SIGNAL
                    ravg1=runningaverage(mean(waveformBrel(23:44,:)),3200);
                    ravg2=runningaverage(mean(waveformBrel(1:22,:)),3200);
                    figure;hold on;plot(ravg1,'b');plot(ravg2,'r')
                    corrcoef(ravg1,ravg2)
                    ravgE=runningaverage(mean(waveformBrel(isESC2,:)),3200);
                    ravgH=runningaverage(mean(waveformBrel(isHIT2,:)),3200);
                    timevals=runningaverage([-2000:1/32:1000],3200);
                    figure;hold on;plot(timevals,ravgE,'b');plot(timevals,ravgH,'r')
                    corrcoef(ravgE,ravgH)
        % SPIKE EST BY THRESHOLD
              clear firingrateB timevals
              for i=1:size(waveformBrel,1)
                  thresh=mean(waveformBrel(i,:))+3*std(waveformBrel(i,:));
                  spikeests=find(waveformBrel(i,:)>thresh);
                  timewinstart=[-2050:20:950];
                  timewinwidth=100;
                  tonset=(fs*2);
                  for j=1:length(timewinstart)
                    firingrateB(i,j)=length(find(spikeests>tonset+timewinstart(j)*32 & spikeests<tonset+timewinstart(j)*32+timewinwidth*32));
                    timevals(j)=timewinstart(j)+timewinwidth/2;
                  end
              end
                figure;hold on;
                plot(timevals,mean(firingrateB(isESC2,:)))
                plot(timevals,mean(firingrateB(isHIT2,:)),'r')

        % SIGNAL AFTER SAUCY
                clear timevals firingrateB
                for i=1:length(spikesBrel)
                  timewinstart=[-2050:20:950];
                    timewinwidth=100;
                    for j=1:length(timewinstart)
                        firingrateB(i,j)=length(find(spikesBrel{i}>timewinstart(j) & spikesBrel{i}<timewinstart(j)+timewinwidth));
                        timevals(j)=timewinstart(j)+timewinwidth/2;
                    end
                end
                figure;hold on;
                plot(timevals,mean(firingrateB(isESC2,:)))
                plot(timevals,mean(firingrateB(isHIT2,:)),'r')

