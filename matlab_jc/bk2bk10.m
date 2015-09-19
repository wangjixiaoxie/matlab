% bk2bk10

% syntax/lesion candidate


% PRELESION
        % started screening 10.04.11
        % decent song
        % A --> BA about 50% of the time
        % A --> C (low stacks) about 50% of the time
        % Template to hit C, does great
    % SYNTAX SHIFT
        % WN on at dawn on 10.06.11 (plan to turn off at dark on sunday)
        % WN off at dark on 10.10.11 - 4 days
        %  60% to 37% 

% 10.18.11 - lesion MMAN & LMAN

% POSTLESION
        % 10.19 - afternoon - began signing again!
        % 10.20 - 92% C, 8% B
        %       - big FF variability reduction
        % eventually recovers to baseline transition probabilities
    % SYNTAX SHIFT
        % 10.24 - wn on for syntax
        % 10.28 - no learning
    % FF SHIFT
        % 11.01 - wn on for FF - hit below 3325Hz (~70%)
                % 2:22pm - hit below 3345Hz
        % 11.02 - wn off at dark
        % 60Hz upshift
        
load DataFF.mat
figure;hold on;
plot(timing3(fvPre),mean(pitchPre(180:230,:)),'.')
plot(timing3(fvWN),mean(pitchWN(180:230,:)),'r.')
plot(timing3(fvPOST),mean(pitchPOST(180:230,:)),'.')
plot(runningaverage(timing3(fvPre),20),runningaverage(mean(pitchPre(210:250,:)),20),'Linewidth',3)
plot(runningaverage(timing3(fvWN),20),runningaverage(mean(pitchWN(210:250,:)),20),'Linewidth',3,'Color','r')
plot(runningaverage(timing3(fvPOST),20),runningaverage(mean(pitchPOST(210:250,:)),20),'Linewidth',3)


nsw=numsortedWN(1100:end);
alls=runningaverage(nsw,12);    
% p=0.022


dirf('*.cbin.not.mat','batchnotes')
fvMANb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvMANc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvMANb=timing3(fvMANb);
tvMANc=timing3(fvMANc);
tvallMAN=timing3([fvMANb fvMANc]);
numsall2=[zeros(1,length(tvMANb)) ones(1,length(tvMANc))];
[b2,ind2]=sort(tvallMAN);
numsortedMAN=numsall2(ind2);
tvsortedMAN=tvallMAN(ind2);

dirf('*.cbin.not.mat','batchnotes')
fvMANWNb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvMANWNc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvMANWNb=timing3(fvMANWNb);
tvMANWNc=timing3(fvMANWNc);
tvallMANWN=timing3([fvMANWNb fvMANWNc]);
numsall2=[zeros(1,length(tvMANWNb)) ones(1,length(tvMANWNc))];
[b2,ind2]=sort(tvallMANWN);
numsortedMANWN=numsall2(ind2);
tvsortedMANWN=tvallMANWN(ind2);

dirf('*.cbin.not.mat','batchnotes')
fvMANPOSTb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvMANPOSTc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvMANPOSTb=timing3(fvMANPOSTb);
tvMANPOSTc=timing3(fvMANPOSTc);
tvallMANPOST=timing3([fvMANPOSTb fvMANPOSTc]);
numsall2=[zeros(1,length(tvMANPOSTb)) ones(1,length(tvMANPOSTc))];
[b2,ind2]=sort(tvallMANPOST);
numsortedMANPOST=numsall2(ind2);
tvsortedMANPOST=tvallMANPOST(ind2);

figure;hold on;
plot(runningaverage([ tvsortedPOSTMAN2 tvsortedPOSTMAN3],40),...
    runningaverage([ numsortedPOSTMAN2 numsortedPOSTMAN3],40),'Color','b')
plot(runningaverage(tvsortedMANWN,40),runningaverage(numsortedMANWN,40),'Color','r')

plot(runningaverage(tvsortedMAN,40),runningaverage(numsortedMAN,40),'Color','b')

figure;hold on;
plot(runningaverage(numsortedMANPOST,40),'Color','r')
plot(runningaverage(numsortedMAN,40),'Color','b')

plot(runningaverage(numsortedMAN,40),'Color','b')


load Data1006.mat

dirf('*.cbin.not.mat','batchnotes')
fvWNb1=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvWNc1=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvWNb=timing3(fvWNb);
tvWNc=timing3(fvWNc);
tvallWN=timing3([fvWNb fvWNc]);
numsall2=[zeros(1,length(tvWNb)) ones(1,length(tvWNc))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);

fvPREb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvPREc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvPREb=timing3(fvPREb);
tvPREc=timing3(fvPREc);
tvallPRE=timing3([fvPREb fvPREc]);
numsall2=[zeros(1,length(tvPREb)) ones(1,length(tvPREc))];
[b2,ind2]=sort(tvallPRE);
numsortedPRE=numsall2(ind2);
tvsortedPRE=tvallPRE(ind2);

fvPOSTb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvPOSTc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvPOSTb=timing3(fvPOSTb);
tvPOSTc=timing3(fvPOSTc);
tvallPOST=timing3([fvPOSTb fvPOSTc]);
numsall2=[zeros(1,length(tvPOSTb)) ones(1,length(tvPOSTc))];
[b2,ind2]=sort(tvallPOST);
numsortedPOST=numsall2(ind2);
tvsortedPOST=tvallPOST(ind2);
% probability of A&C as a function of time (i.e. not B)
figure;hold on;
plot(runningaverage(tvsortedPRE,80),runningaverage(numsortedPRE,80),'Color','b')
plot(runningaverage(tvsortedWN,80),runningaverage(numsortedWN,80),'Color','r')
plot(runningaverage(tvsortedPOST,50),runningaverage(numsortedPOST,50),'Color','b')

fvPOSTMAN2b=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvPOSTMAN2c=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvPOSTMAN2b=timing3(fvPOSTMAN2b);
tvPOSTMAN2c=timing3(fvPOSTMAN2c);
tvallPOSTMAN2=timing3([fvPOSTMAN2b fvPOSTMAN2c]);
numsall2=[zeros(1,length(tvPOSTMAN2b)) ones(1,length(tvPOSTMAN2c))];
[b2,ind2]=sort(tvallPOSTMAN2);
numsortedPOSTMAN2=numsall2(ind2);
tvsortedPOSTMAN2=tvallPOSTMAN2(ind2);



%%%%%%%%%%%%%
%%%%%%%%%%%%%

% postMANlesion

% 1019_postMANlesion
                songstarts=find(labels=='/');
                songstats=[];
                count=0;
                for i=1:length(songstarts)-1
                    bps=find(labels(songstarts(i):songstarts(i+1))=='a');
                    % if there is a branch point
                    if ~isempty(bps)
                        songstats=[songstats -1]; % song begins
                        for j=1:length(bps)
                            if labels(songstarts(i)+bps(j))=='c'
                                nxt=1;
                            else
                                nxt=0;
                            end
                            songstats=[songstats nxt];
                        end
                    end
                end

        % Actual number of runs
            from0to1=find(songstats(1:end-1)==0 & songstats(2:end)==1);
            from1to0=find(songstats(1:end-1)==1 & songstats(2:end)==0);
            actualnumruns=length(from0to1)+length(from1to0);
            % Actual number is 285 - p<0.01 - see below - very high number of runs
            % suggests propensity for alternation

        % Resampled number of runs
            songstats3=songstats(songstats~=-1);
            pnull=mean(songstats3);
            randomdraws=rand(1,1e6);
            randbindraws=randomdraws<pnull;
            count=0;
            clear numruns
            for z=1:1000
                % make a mock data set
                clear mockss
                for i=1:length(songstats)
                    if songstats(i)==1 | songstats(i)==0
                        count=count+1;
                        mockss(i)=randbindraws(count);
                    else
                        mockss(i)=-1;
                    end
                end
                % calculate the number of runs
                from0to1=find(mockss(1:end-1)==0 & mockss(2:end)==1);
                from1to0=find(mockss(1:end-1)==1 & mockss(2:end)==0);
                numruns(z)=length(from0to1)+length(from1to0);
            end
        length(find(numruns>actualnumruns-1))/length(numruns)
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% 1004_tmptest
        labels=notestats
                songstarts=find(labels=='/');
                songstats=[];
                count=0;
                for i=1:length(songstarts)-1
                    bps=find(labels(songstarts(i):songstarts(i+1))=='a');
                    % if there is a branch point
                    if ~isempty(bps)
                        songstats=[songstats -1]; % song begins
                        for j=1:length(bps)
                            if labels(songstarts(i)+bps(j))=='c'
                                nxt=1;
                            else
                                nxt=0;
                            end
                            songstats=[songstats nxt];
                        end
                    end
                end
                
    from0to1=find(songstats(1:end-1)==0 & songstats(2:end)==1);
    from1to0=find(songstats(1:end-1)==1 & songstats(2:end)==0);
    actualnumruns=length(from0to1)+length(from1to0);
    % Actual number is 334 - p<0.0001 - see below - very high number of runs
    % suggests propensity for alternation

% Resampled number of runs
            songstats3=songstats(songstats~=-1);
            pnull=mean(songstats3);
            randomdraws=rand(1,1e6);
            randbindraws=randomdraws<pnull;
            count=0;
            clear numruns
            for z=1:1000
                % make a mock data set
                clear mockss
                for i=1:length(songstats)
                    if songstats(i)==1 | songstats(i)==0
                        count=count+1;
                        mockss(i)=randbindraws(count);
                    else
                        mockss(i)=-1;
                    end
                end
                % calculate the number of runs
                from0to1=find(mockss(1:end-1)==0 & mockss(2:end)==1);
                from1to0=find(mockss(1:end-1)==1 & mockss(2:end)==0);
                numruns(z)=length(from0to1)+length(from1to0);
            end
        length(find(numruns>actualnumruns-1))/length(numruns)
    