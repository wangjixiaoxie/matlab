% Repeat analysis
% [timevals,distnt] = jcrepeatdist('a','batchnotes');
% cardinal1/Evren_Repeats/README - read this


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% June 20, 2011
% More detailed behavioral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clear all
        load /cardinal9/Repeat_Analysis/Exp_062011.mat
% For each experiment

%
for i=1:length(Exp)
    inder(i)=Exp(i).postlesion==0 & Exp(i).hitabove==1;
end
ind1=find(inder);
% no lesion, learning to shift down
% For each experiment, calculate ptrans, repeat lengths, and hit rates
for i=ind1;
    MINI=Exp(i).MIN;
    MINIless=MINI-1;
    clear lastptrans replong hitrate day
    count=0;
    clear semptrans semreplong semhitrate
    for j=1:3 % Pre, Dur, Post
            if Exp(i).JC==0
                unidays=find(diff(Exp(i).time{j})>0.3); % Evren time
            else
                unidays=find(diff(Exp(i).time{j})>8); % JC time
            end
            unistart=[1 unidays+1];
            uniend=[unidays length(Exp(i).time{j})];
            if i==6 & j==3 % Missing data for this experiment (Evren)
                lastptrans(j).data(1)=0;
                day(j).data(1)=0;
                replong(j).data(1)=0;
                hitrate(j).data(1)=0;
                semptrans(j).data(1)=0;
                semreplong(j).data(1)=0;
                semhitrate(j).data(1)=0;
            else

            for ii=1:length(unistart) % For each day
                rpdata=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                count=count+1;
                lastptrans(j).data(ii)=length(find(rpdata>=MINI))/length(find(rpdata>=MINIless));
                day(j).data(ii)=count;
                replong(j).data(ii)=mean(rpdata);
                hitrate(j).data(ii)=sum(rpdata(find(rpdata>MINIless))-MINIless*ones(1,length(find(rpdata>MINIless))))/sum(rpdata);
                % Create surrogate data sets using resampling
                clear ptransBS replongBS hitrateBS
                for bs=1:1000
                    scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                    scrambleRPL=rpdata(scrambleIND);
                    ptransBS(bs)=length(find(scrambleRPL>=MINI))/length(find(scrambleRPL>=MINIless));
                    replongBS(bs)=mean(scrambleRPL);
                    hitrateBS(bs)=sum(scrambleRPL(find(scrambleRPL>MINIless))-MINIless*ones(1,length(find(scrambleRPL>MINIless))))/sum(scrambleRPL);
                end
                % The standard deviation of these surrogate data sets is equivalent
                % to the standard error of the mean for the statistic of
                % interest
                semptrans(j).data(ii)=std(ptransBS);
                semreplong(j).data(ii)=std(replongBS);
                semhitrate(j).data(ii)=std(hitrateBS);
            end
        end
    end
    [probstay,numstay,semstay]=repeatsTOptrans2(Exp,i);
    % Resample these probabilities

    Exp(i).probstay=probstay;
    Exp(i).numstay=numstay;
    Exp(i).semstay=semstay;
    Exp(i).lastptrans=lastptrans;
    Exp(i).day=day;
    Exp(i).replong=replong;
    Exp(i).hitrate=hitrate;
    Exp(i).semptrans=semptrans;
    Exp(i).semreplong=semreplong;
    Exp(i).semhitrate=semhitrate;
end

% Transition probability reduction
        % Point 2 is final day with WN
        % All 
        figure;hold on;
        for i=ind1
            indNN=find(~isnan(Exp(i).probstay(1).rn) & ~isnan(Exp(i).probstay(2).rn) & Exp(i).numstay(1).rn >5 & Exp(i).numstay(2).rn>5);
            plot(Exp(i).probstay(1).rn(indNN),Exp(i).probstay(2).rn(indNN),'b.','Markersize',15)
            plot([Exp(i).probstay(1).rn(indNN);Exp(i).probstay(1).rn(indNN)],[Exp(i).probstay(2).rn(indNN)-Exp(i).semstay(2).rn(indNN);Exp(i).probstay(2).rn(indNN)+Exp(i).semstay(2).rn(indNN)],'b-')
            plot([Exp(i).probstay(1).rn(indNN)-Exp(i).semstay(1).rn(indNN);Exp(i).probstay(1).rn(indNN)+Exp(i).semstay(1).rn(indNN)],[Exp(i).probstay(2).rn(indNN);Exp(i).probstay(2).rn(indNN)],'b-')
        end
        % Transition immediately before targeting
        for i=ind1
            indNN=Exp(i).MIN-1;
            plot(Exp(i).probstay(1).rn(indNN),Exp(i).probstay(2).rn(indNN),'r.','Markersize',25)
            plot([Exp(i).probstay(1).rn(indNN);Exp(i).probstay(1).rn(indNN)],[Exp(i).probstay(2).rn(indNN)-Exp(i).semstay(2).rn(indNN);Exp(i).probstay(2).rn(indNN)+Exp(i).semstay(2).rn(indNN)],'r-')
            plot([Exp(i).probstay(1).rn(indNN)-Exp(i).semstay(1).rn(indNN);Exp(i).probstay(1).rn(indNN)+Exp(i).semstay(1).rn(indNN)],[Exp(i).probstay(2).rn(indNN);Exp(i).probstay(2).rn(indNN)],'r-')
        end
        plot([0 1],[0 1],'k')
% Hit rate reduction
figure;hold on;
        for i=ind1
            plot(mean(Exp(i).hitrate(1).data),Exp(i).hitrate(2).data(end),'b.','Markersize',15)
            mnpre=[mean(Exp(i).hitrate(1).data);mean(Exp(i).hitrate(1).data)]
            sempost=[mean(Exp(i).hitrate(2).data)-mean(Exp(i).semhitrate(2).data);mean(Exp(i).hitrate(2).data)+mean(Exp(i).semhitrate(2).data)];
            mnpost=[mean(Exp(i).hitrate(2).data);mean(Exp(i).hitrate(2).data)];
            sempre=[mean(Exp(i).hitrate(1).data)-mean(Exp(i).semhitrate(1).data);mean(Exp(i).hitrate(1).data)+mean(Exp(i).semhitrate(1).data)];
            plot(mnpre,sempost,'b-')
            plot(sempre,mnpost,'b-')
        end
        plot([0 1],[0 1],'k')



    figure;hold on;
    subplot(311);hold on;
    for i=1:3
        if i==2; clrvl='r'; else clrvl='k';end
        plot(thisday(i).data,thisreplong(i).data,'Color',clrvl,'Marker','.','Markersize',15)
        plot([thisday(i).data;thisday(i).data],[thisreplong(i).data+semrplong(i).data;thisreplong(i).data-semrplong(i).data],'Color',clrvl)
    end
    subplot(312);hold on;
    for i=1:3
        if i==2; clrvl='r'; else clrvl='k';end
        plot(thisday(i).data,thisptrans(i).data,'Color',clrvl,'Marker','.','Markersize',15)
        plot([thisday(i).data;thisday(i).data],[thisptrans(i).data+semptrans(i).data;thisptrans(i).data-semptrans(i).data],'Color',clrvl)
    end
    subplot(313);hold on;
    for i=1:3
        if i==2; clrvl='r'; else clrvl='k';end
        plot(thisday(i).data,thishitrate(i).data,'Color',clrvl,'Marker','.','Markersize',15)
        plot([thisday(i).data;thisday(i).data],[thishitrate(i).data+semhitrate(i).data;thishitrate(i).data-semhitrate(i).data],'Color',clrvl)
    end

end
    
    

    
        




% 1. Repeat number change

% 2. Change in hit rate

% 3. Reduction in targeted p(trans)

% 4. Reduction in p(trans) relative to baseline p(trans)


        
        %%%%%%%%%%
        % REPEAT NUMBER - summarize change before and after lesion as V plots (pre-dur-post) to show reduction during WN
        %%%%%%%%%%
        [alls,alls2]=plotRepeatLearningSummary1(Exp)
        % first plot includes all three prelesion exps for w93pk62
        % second plot only includes the one exp that I did (ignores Evren's)
        % alls2 - three columns for pre,dur,post (each is the last day of the period) 
        %       - n rows corresponding to each entry in Exp (i.e. each experiment)

        %%%%%%%%%%%%%%%
        % Calculate transition probabilities (probstay)
        % Note that numstay corresponds to sample sizes
        %%%%%%%%%%%%%%%
        [probstay,numstay]=repeatsTOptrans(Exp);
        % returns the following:
        % probstay(experiment,[pre=1,dur=2,post=3]).rn -- where rn(i) is the
            % probability of the ith transition, i.e. the transition from the ith to
            % the i+1th rendition of the repeating syllable
        % e.g. probstay(11,2).rn(1) is the probability of the transition from the
            % 1st to the 2nd rendition of the syllable in the 11th experiment during
            % the last day of WN 
        %%%%%%%%%%%%%%%%%%%
        % Plot all learning data pre and post-lesion as both p(transition) and repeat numbers
        %%%%%%%%%%%%%%%%%
        indPost=[1 2 14 15];
        plotRepeatLearningSummary2(Exp,probstay,numstay,alls2,indPost)
        
  % Give data to Tim
  clear TWdata
        for i=1:15
            TWdata(i).ptrans{1}=probstay(i,1).rn;TWdata(i).ptrans{2}=probstay(i,2).rn;TWdata(i).ptrans{3}=probstay(i,3).rn;
            TWdata(i).sampsize{1}=numstay(i,1).rn;TWdata(i).sampsize{2}=numstay(i,2).rn;TWdata(i).sampsize{3}=numstay(i,3).rn; 
            TWdata(i).rplength=alls2(i,:);
            TWdata(i).MINptrans_hit=Exp(i).MIN-1;
            TWdata(i).bird=Exp(i).bird;
            TWdata(i).postlesionkeepind=[1 2 14 15];
            TWdata(i).prelesionkeepind=[4 6 10 12];
        end

sum(distntWN(find(distntWN>4))-4*ones(1,length(find(distntWN>4))))/sum(distntWN)
%%%%%%%
[probstay,numstay]=repeatsTOptrans(Exp);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% June 3, 2011 
% - stationarity over a month example figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% April 19, 2011 
% - stationarity over a months
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /cardinal9/SyntaxBirds/stationaryMONTH.mat
for j=1:length(Rep)
    for i=1:max([Rep(j).rplengthPre Rep(j).rplengthPost])
        Rep(j).ptransPre(i)=length(find(Rep(j).rplengthPre>=i+1))/length(find(Rep(j).rplengthPre>=i));
        Rep(j).num_in_ptransPre(i)=length(find(Rep(j).rplengthPre>=i+1));
        Rep(j).ptransPost(i)=length(find(Rep(j).rplengthPost>=i+1))/length(find(Rep(j).rplengthPost>=i));
        Rep(j).num_in_ptransPost(i)=length(find(Rep(j).rplengthPost>=i+1));       
    end
end
figure;hold on;
min_n=15;
plot(ptransPre(find(replongPre>min_n & replongPost>min_n)),ptransPost(find(replongPre>min_n & replongPost>min_n)),'.')
plot([0 1],[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% April 12, 2011 
% - figures for lab meeting - transition probabilities during learning %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear all
        load /cardinal9/SyntaxLearningFigs/Repeats5.mat
        %%%%%%%%%%
        % REPEAT NUMBER - summarize change before and after lesion as V plots (pre-dur-post) to show reduction during WN
        %%%%%%%%%%
        [alls,alls2]=plotRepeatLearningSummary1(Exp)
        % first plot includes all three prelesion exps for w93pk62
        % second plot only includes the one exp that I did (ignores Evren's)
        % alls2 - three columns for pre,dur,post (each is the last day of the period) 
        %       - n rows corresponding to each entry in Exp (i.e. each experiment)

        %%%%%%%%%%%%%%%
        % Calculate transition probabilities (probstay)
        % Note that numstay corresponds to sample sizes
        %%%%%%%%%%%%%%%
        [probstay,numstay]=repeatsTOptrans(Exp);
        % returns the following:
        % probstay(experiment,[pre=1,dur=2,post=3]).rn -- where rn(i) is the
            % probability of the ith transition, i.e. the transition from the ith to
            % the i+1th rendition of the repeating syllable
        % e.g. probstay(11,2).rn(1) is the probability of the transition from the
            % 1st to the 2nd rendition of the syllable in the 11th experiment during
            % the last day of WN 
        %%%%%%%%%%%%%%%%%%%
        % Plot all learning data pre and post-lesion as both p(transition) and repeat numbers
        %%%%%%%%%%%%%%%%%
        indPost=[1 2 14 15];
        plotRepeatLearningSummary2(Exp,probstay,numstay,alls2,indPost)
        
  % Give data to Tim
  clear TWdata
        for i=1:15
            TWdata(i).ptrans{1}=probstay(i,1).rn;TWdata(i).ptrans{2}=probstay(i,2).rn;TWdata(i).ptrans{3}=probstay(i,3).rn;
            TWdata(i).sampsize{1}=numstay(i,1).rn;TWdata(i).sampsize{2}=numstay(i,2).rn;TWdata(i).sampsize{3}=numstay(i,3).rn; 
            TWdata(i).rplength=alls2(i,:);
            TWdata(i).MINptrans_hit=Exp(i).MIN-1;
            TWdata(i).bird=Exp(i).bird;
            TWdata(i).postlesionkeepind=[1 2 14 15];
            TWdata(i).prelesionkeepind=[4 6 10 12];
        end
%%%%%%%%%%%%%%%%%%%%
%%% April 7, 2011 %%
%%%%%%%%%%%%%%%%%%%%
% Figures for Tim - repeat learning example (pre-lesion)
%%%%%%%%%%%%%%%%%%%%
    load /cardinal9/SyntaxBirds/LearningExample.mat

    % from /cardinal1/Evren_Repeats/w93pk62/ampon/batch.notcatch
    [avZ,t,f]=get_avn('batch.notcatch','z',0.3,3.5,'-','','obs0');
    figure
    imagesc(t,f,log(avZ));syn;
    % 
    figure;hold on;
    plot(runningaverage(timevalsDUR,20),runningaverage(distntDUR,20),'.')
    plot(runningaverage(timevalsPRE,20),runningaverage(distntPRE,20),'k.')
    plot(runningaverage(timevalsPOST,20),runningaverage(distntPOST,20),'k.')
 
            time{1}=timevalsPRE;
            time{2}=timevalsDUR;
            time{3}=[timevalsPOST timevals1];
            rplength{1}=distntPRE;
            rplength{2}=distntDUR;
            rplength{3}=[distntPOST distnt1];
            MINI=3;
            for j=1:3
                unidays=find(diff(time{j})>8); % JC time
                unistart=[1 unidays+1];
                uniend=[unidays length(time{j})];
                for k=MINI-1
                    for i=1:length(unistart)
                        ptrans(i,j)=length(find(rplength{j}(unistart(i):uniend(i))>=k+1))/length(find(rplength{j}(unistart(i):uniend(i))>=k));
                        replong(i,j)=mean(rplength{j}(unistart(i):uniend(i)));
                    end
                end
            end
            figure;hold on;
            subplot(212);hold on;
            plot([1:9],[ptrans(1,1);ptrans(2,1);ptrans(1:4,2);ptrans(1:3,3)],'.','Markersize',15)
            plot([1:9],[ptrans(1,1);ptrans(2,1);ptrans(1:4,2);ptrans(1:3,3)],'-')
            plot([1:9],1-[ptrans(1,1);ptrans(2,1);ptrans(1:4,2);ptrans(1:3,3)],'r.','Markersize',15)
            plot([1:9],1-[ptrans(1,1);ptrans(2,1);ptrans(1:4,2);ptrans(1:3,3)],'r-')
            plot([0 10],[mean(ptrans(1:2,1)) mean(ptrans(1:2,1))])
            plot([0 10],1-[mean(ptrans(1:2,1)) mean(ptrans(1:2,1))])
            ylim([0 1])
            subplot(211);hold on;
            plot([1:9],[replong(1,1);replong(2,1);replong(1:4,2);replong(1:3,3)],'.','Markersize',15)
            plot([1:9],[replong(1,1);replong(2,1);replong(1:4,2);replong(1:3,3)],'-')
            plot([0 10],[mean(replong(1:2,1)) mean(replong(1:2,1))])
            ylim([2 5])

    
        
%%%%%%%%%%%%%%%%%%%%
%%% April 5, 2011 %%
%%%%%%%%%%%%%%%%%%%%
% Acute effects of lesions and inactivations on repeat number
%%%%%%%%%%%%%%%%%%%%
        clear all
        load /cardinal3/SyntaxBirds/LesionsANDinas0404.mat
        Experiment=[Exps1 Exps2 Exps3 Exps4];

        % Plot Summary data and calculate variables for further analysis
        % plot colors: black - Evlesions, green - Ourlesions, blue - muscimol in LMAN, red - APV in RA
        [mnpre,sdpre,mnpost,sdpost,indlesion,indours,indMU,indAPV,notetype]=plotAcuteLesionEffectsonRepeats(Experiment);

        % How to explain the variation? Relationship to baseline repeat count? Not much
            figure;hold on;
            msize=15;
            subplot(321);hold on;
            plot(mnpre(indlesion),mnpost(indlesion)./mnpre(indlesion),'k.','Markersize',msize)
            plot(mnpre(indAPV),mnpost(indAPV)./mnpre(indAPV),'r.','Markersize',msize)
            subplot(322);hold on;
            plot(mnpre(indlesion),mnpost(indlesion)-mnpre(indlesion),'k.','Markersize',msize)
            plot(mnpre(indAPV),mnpost(indAPV)-mnpre(indAPV),'r.','Markersize',msize)
        % How to explain the variation? Relationship to baseline repeat CV? No
            subplot(323);hold on;
            plot(sdpre(indAPV)./mnpre(indAPV),mnpost(indAPV)./mnpre(indAPV),'r.','Markersize',msize)
            plot(sdpre(indlesion)./mnpre(indlesion),mnpost(indlesion)./mnpre(indlesion),'k.','Markersize',msize)   
            subplot(324);hold on;
            plot(sdpre(indAPV),mnpost(indAPV)./mnpre(indAPV),'r.','Markersize',msize)
            plot(sdpre(indlesion),mnpost(indlesion)./mnpre(indlesion),'k.','Markersize',msize)   
        % How to explain the variation? note type? not really
            figure;hold on;msize=25;
            % none % plot(1,mnpost(find(notetype(indAPV)==1))./mnpre(find(notetype(indAPV)==1)),'k.','Markersize',msize)
            plot(1,mnpost(find(notetype(indAPV)==2))./mnpre(find(notetype(indAPV)==2)),'b.','Markersize',msize)
            plot(1,mnpost(find(notetype(indAPV)==3))./mnpre(find(notetype(indAPV)==3)),'r.','Markersize',msize)    
            plot(1,mnpost(find(notetype(indAPV)==4))./mnpre(find(notetype(indAPV)==4)),'g.','Markersize',msize)    
            msize=25;
            plot(2,mnpost(find(notetype(indlesion)==1))./mnpre(find(notetype(indlesion)==1)),'k.','Markersize',msize)
            plot(2,mnpost(find(notetype(indlesion)==2))./mnpre(find(notetype(indlesion)==2)),'b.','Markersize',msize)
            plot(2,mnpost(find(notetype(indlesion)==3))./mnpre(find(notetype(indlesion)==3)),'r.','Markersize',msize)    
            % none %plot(2,mnpost(find(notetype(indlesion)==4))./mnpre(find(notetype(indlesion)==4)),'g.','Markersize',msize)        
            xlim([0.5 2.5])

        % How to explain the variation? FF variability reduction? NO
            indpitch=find(ptmeasureind(:,1)>0);
            for k=1:length(indpitch)
                j=indpitch(k);
                sdratio(j)=mean(jcstd(Experiment(j).pitchMU(ptmeasureind(j,1):ptmeasureind(j,2),:)')./jcstd(Experiment(j).pitchACpre(ptmeasureind(j,1):ptmeasureind(j,2),:)'));
            end
            figure;plot(1-sdratio(indpitch),mnpost(indpitch)./mnpre(indpitch),'.','Markersize',msize) % no relationship
 
         % Same bird, different effects
        figure;hold on;
        subplot(321);hold on;
        [b,a]=hist(Experiment(24).distnt{1},[1:1:12]);
        [c,d]=hist(Experiment(24).distnt{2},[1:1:12]);
        plot(b/sum(b),'b');plot(c/sum(c),'r')
        subplot(322);hold on;
        [b,a]=hist(Experiment(25).distnt{1},[1:1:12]);
        [c,d]=hist(Experiment(25).distnt{2},[1:1:12]);
        plot(b/sum(b),'b');plot(c/sum(c),'r')
        
      load /cardinal3/SyntaxBirds/bk35pitch.mat
subplot(323);hold on;plot(pitchpre24(:,1:30),'b');plot(pitchpost24(:,1:30),'r')
xlim([200 450]);ylim([2700 3600])
subplot(324);hold on;plot(pitchpre25(:,1:30),'b');plot(pitchpost25(:,1:30),'r')
xlim([150 280]);ylim([1000 2000])
subplot(325);hold on;plot(std(pitchpre24'),'b');plot(std(pitchpost24'),'r')
xlim([200 450]);ylim([0 110])
subplot(326);hold on;plot(std(pitchpre25'),'b');plot(std(pitchpost25'),'r')
xlim([150 280]);ylim([0 130])
%%%%%%%%%%%%%%%%%%%%
%%% March 11, 2011 %%
%%%%%%%%%%%%%%%%%%%%        
% STATIONARITY of baseline data 
% Looked at Evren's data & some of my data from other folders (e.g. bk5bk25
% Generally seems easy to find 4-5 days in a row.
    % First day of recording may produce outliers...?
                load /cardinal9/SyntaxBirds/stationaryFinal.mat
                % n=9 in n=6 birds
                figure;hold on;
                count=0;
                for i=1:length(stationary) % for each bird
                    for j=1:length(stationary(i).timevals) % for each repeated syllable
                        count=count+1;
                        subplot(4,3,count)
                        plot(runningaverage(stationary(i).timevals{j}-min(stationary(i).timevals{j}),2)/24,...
                            runningaverage(stationary(i).distnt{j},2),'.')
                        ylim([0 14])
                    end
                end
                count=0;
                clear probstay numstay avgrplength
                for i=1:length(stationary) % for each bird
                    for j=1:length(stationary(i).timevals) % for each repeated syllable
                        count=count+1;
                        mntime=min(stationary(i).timevals{j});
                        for k=1:5
                            ind1=find(stationary(i).timevals{j}>mntime+(k-1)*24 & stationary(i).timevals{j}<mntime+k*24);
                            avgrplength(count,k)=mean(stationary(i).distnt{j}(ind1));
                            % ptrans
                            for m=2:20
                                probstay(count).day(m-1,k)=length(find(stationary(i).distnt{j}(ind1)>m))/length(find(stationary(i).distnt{j}(ind1)>m-1));
                                numstay(count).day(m-1,k)=length(find(stationary(i).distnt{j}(ind1)>m));
                            end
                        end
                    end
                end
              
              %%%%%%%%%%%
              %%%%%%%%%%%%
              % Example histograms
              figure;hold on;
                for i=6
                    for j=1
                        mntime=min(stationary(i).timevals{j});
                        for k=1:5
                            ind1=find(stationary(i).timevals{j}>mntime+(k-1)*24 & stationary(i).timevals{j}<mntime+k*24);
                            avgexample(k)=mean(stationary(i).distnt{j}(ind1));
                            for kk=1:500 % resampling to obtain CIs
                                nrsamp=ind1(ceil(length(ind1)*rand(1,length(ind1))));
                                msmn(kk)=mean(stationary(i).distnt{j}(nrsamp));
                            end
                            prchigh(k)=prctile(msmn,97.5); % 95prct
                            prclow(k)=prctile(msmn,2.5);
                        end
                    end
                end
              %%%%%%%%%%%%%%%
              %%%%%%%%%%%%%%%
              % EXAMPLE SPECTROGRAM - bk34bk33/1203_screen - tafsim for 'z'
              % EXAMPLE HISTOGRAMS
              clear d
              for i=1:5
                  ind1=find(stationary(6).timevals{1}>mntime+(i-1)*24 & stationary(6).timevals{1}<mntime+i*24);
                  [b,a]=hist(stationary(6).distnt{1}(ind1),[0:1:15]);
                  [c,d(i,:)]=stairs(a,b/sum(b))
              end
              figure;plot(c,d,'Linewidth',2)
              % EXAMPLE FIGURE - also bk34bk33
                  figure;hold on;
                  for i=1:5
                      plot([i i],[prclow(i) prchigh(i)])
                      plot(i,avgexample(i),'.','Markersize',20)
                  end
                  plot(avgexample)
                  plot([1 5],[mean(avgexample) mean(avgexample)])
                  xlim([0 5.2]);ylim([0 10])
             %%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%
             %%% All data - plot day one vs day four, n=9
              figure;hold on
              nn=9
                for i=1:nn
                    plot(avgrplength(i,1),avgrplength(i,4),'ko','Markersize',5)  
                end    
                plot([3:1:15],[3:1:15],'k')
                ylim([3 13]);xlim([3 13])
             %%%% DO IT FOR PTRANS 
             figure;hold on;
             nn=9
             for i=1:nn
                 for j=1:length(probstay(i).day(:,1))
                     % - note that I only keep the ones with a decently large sample size
                     if mean(numstay(i).day(j,1:4))>10
                        plot(probstay(i).day(j,1),probstay(i).day(j,4),'ko','Markersize',5)
                     end
                 end
             end
             for i=1:nn
                 for j=1:length(probstay(i).day(:,1))
                     % - note that I only keep the ones with a decently large sample size
                     if mean(numstay(i).day(j,1:4))>25
                        plot(probstay(i).day(j,1),probstay(i).day(j,4),'ro','Markersize',5)
                     end
                 end
             end
             plot([0 1],[0 1],'k')
             ylim([0 1]);xlim([0 1])

                
                
                
             %%% All data - means, n=9
             nn=9
             figure;hold on;
             for i=1:nn
                 residrplength(i,1:4)=avgrplength(i,1:4)-mean(avgrplength(i,1:4));
             end
             figure;hold on;
             bar([1:4],mean(abs(residrplength)),'Facecolor','w','Edgecolor','k')
             for i=1:4
                 plot([i i],[mean(abs(residrplength(:,i)))-std(abs(residrplength(:,i)))/sqrt(nn)...
                     mean(abs(residrplength(:,i)))+std(abs(residrplength(:,i)))/sqrt(nn)],'k')
             end
                % Other representations
                figure;hold on; % n=5 birds, n=9 repeats
                subplot(221);hold on;
                for i=1:9
                    plot(avgrplength(i,:)/mean(avgrplength(i,:)))
                end
                ylim([0.8 1.2])
                subplot(222);hold on;
                for i=1:9
                    plot(avgrplength(i,:)-mean(avgrplength(i,:)))
                end    
                ylim([-1 1])
                subplot(223);hold on;
                for i=1:9
                    plot(avgrplength(i,:))
                end    
                subplot(224);hold on;
                for i=1:11
                    plot(avgrplength(i,1),avgrplength(i,2),'k.','Markersize',15)
                    plot(avgrplength(i,1),avgrplength(i,3),'b.','Markersize',15)   
                    plot(avgrplength(i,1),avgrplength(i,4),'r.','Markersize',15)  
                end    
                plot([3:1:15],[3:1:15
                ylim([3 13]);xlim([3 13])

% ACUTE effects of LESION
    % get bk1bk29 and bk67bk44 pre-data from TW
    % lots of Evren data on back-up drive
    % See notes in black marble notebook
    load /cardinal3/LesionsANDinas.mat
    
% ACUTE effects of APV    
    
    
    
%%%%%%%%%%%%%%%%%%%%
%%% March 9, 2011 %%
%%%%%%%%%%%%%%%%%%%%        
Exp=Exp1;
inddown=[1 2 4 6 7 9 10 12 14 15];% downshifts  
clear transition
        for i=inddown
                    for j=1:3 % pre,2nd day of wn
                        if i==6 & j==3
                            break
                        end
                        if Exp(i).JC==1
                            unidays=find(diff(Exp(i).time{j})>8); % JC time
                            unistart=[1 unidays+1];
                            uniend=[unidays length(Exp(i).time{j})];
                        else
                            unistart=[];
                            uniend=[];
                            unidays=unique(floor(Exp(i).time{j})); % Ev time
                            for kk=1:length(unidays)
                                unistart(kk)=min(find(floor(Exp(i).time{j})==unidays(kk)));
                                uniend(kk)=max(find(floor(Exp(i).time{j})==unidays(kk)));
                            end
                        end
                        clear remain
                        if j==2 % second half of second day
                            for k=1:Exp(i).MIN+5
                                length2=floor(uniend(end)-unistart(end)/2)-unistart(end);
                                remain(k)=length(find(Exp(i).rplength{j}(unistart(end)+length2:uniend(end))>=k));
                            end
                        else % entire pre day
                            for k=1:Exp(i).MIN+5
                                remain(k)=length(find(Exp(i).rplength{j}(unistart(end):uniend(end))>=k));
                            end
                        end
                        for k=1:length(remain)-1
                            % given that you've done the kth repeat, what
                            % is chance you'll perseverate 
                            transition(i).pstay(j,k)=remain(k+1)/remain(k);
                            transition(i).pleave(j,k)=1-transition(i).pstay(j,k);
                            transition(i).size(j,k)=remain(k+1);
                            
                        end
                    end
        end
        for i=1:length(Exp)
            transition(i).MIN=Exp(i).MIN-1;
        end
        % "Residuals" - i.e. changes
            for i=inddown
                for j=2:3
                    if i==6 & j==3
                        break
                    end
                    for k=1:length(transition(i).pstay(j,:))
                        ind=10-transition(i).MIN+k;
                        transition(i).pstayresids(j,ind)=transition(i).pstay(j,k)-transition(i).pstay(1,k);
                        transition(i).goodind(j,ind)=transition(i).size(j,k)>2;
                    end
                end
            end
            figure;hold on;
            for i=inddown
                plot(find(transition(i).goodind(2,:))-10+1,transition(i).pstayresids(2,transition(i).goodind(2,:))','.','Markersize',15,'color','b')
                if size(transition(i).pstayresids,1)>2 % RECOVERY
                    %plot(find(transition(i).goodind(3,:))-10+1,transition(i).pstayresids(3,transition(i).goodind(3,:))','.','Markersize',15,'color','r')
                end
                if Exp(i).postlesion
                plot(find(transition(i).goodind(2,:))-10+1,transition(i).pstayresids(2,transition(i).goodind(2,:))','.','Markersize',15,'color','r')
                end

            end
            plot([-3.5 3.5],[0 0],'k')
            xlim([-3.5 3.5])
   transitionTW=transition([1 2 4 6 7 9 10 12 14]);
     
        
  % PLOT - note that the
        figure;hold on;
        for i=1:length(inddown)
            subplot(4,4,i);hold on;
            indexthis=inddown(i);
            plot(transition(indexthis).pstay')
            plot(transition(indexthis).pstay','.')            
            plot([Exp(indexthis).MIN-1 Exp(indexthis).MIN-1],[0 1],'k')
%             for j=1:3
%             transition(i).delt(j)=transition(i).pstay
%             end
        end
            
%        
        
%%%%%%%%%%%%%%%%%%%%
%%% March 1, 2011 %%
%%%%%%%%%%%%%%%%%%%%

% Note - FOR NUMBER 6, there is no post data
% Exp(i).averages(1) ---> average repeat length for day prior to wn
% Exp(i).averages(2) ---> average repeat length for second half of all
                    % catch trials sung on last day of wn (typically day 2)
% Exp(i).averages(3) ---> average repeat length for day n+2 (where day n is
                    % last day of wn)
                    
        for i=[1:5 7:15]%length(Exp) % #6 - no post data
                prenum=Exp(i).pre;
                    mintime=min(Exp(i).time{prenum});
                    for j=1:length(Exp(i).time)
                        if Exp(i).JC==1
                            unidays=find(diff(Exp(i).time{j})>8); % JC time
                            unistart=[1 unidays+1];
                            uniend=[unidays length(Exp(i).time{j})];
                        else
                            unistart=[];
                            uniend=[];
                            unidays=unique(floor(Exp(i).time{j})); % Ev time
                            for kk=1:length(unidays)
                                unistart(kk)=min(find(floor(Exp(i).time{j})==unidays(kk)));
                                uniend(kk)=max(find(floor(Exp(i).time{j})==unidays(kk)));
                            end
                        end
                        if j==1 % day before wn
                            Exp(i).averages(j)=mean(Exp(i).rplength{j}(unistart(end):uniend(end)));
                        else if j==2 % second half of last day with wn
                                length2=floor(uniend(end)-unistart(end)/2)-unistart(end);
                            Exp(i).averages(j)=mean(Exp(i).rplength{j}(unistart(end)+length2:uniend(end)));
                        else if j==3 % day n+2 or more, where day n is last day with wn
                            Exp(i).averages(j)=mean(Exp(i).rplength{j}(unistart(end):uniend(end)));
                        end
                        end
                        end
                    end
        end
%
%%%% 
% 
figure;hold on;
for i=[1:length(Exp)] % only one up and down point for w93pk62
    if Exp(i).postlesion==0
    if i~=6
        plot([1 2 3],[Exp(i).averages])
    else
        plot([1 2],[Exp(i).averages])
    end
    else
        plot([1 2 3],[Exp(i).averages],'r')
    end
end



%%%%%%%%%%5
%%%%%%%%%%%%
%%%%%%%%%%%% 



        clear all
        load /cardinal3/SyntaxLearningFigs/Repeats4.mat
        Exp=Exp1;
        figure;hold on;
        lenham=11;
        h=hamming(lenham);
        h=h*(lenham/sum(h));

        for i=1:length(Exp) %- 2 is bk1bk29 - better example
            if i<4
                lenham=51;
                h=hamming(lenham);
                h=h*(lenham/sum(h));
            end
            if i>3
                lenham=11;
                h=hamming(lenham);
                h=h*(lenham/sum(h));
            end
            subplot(4,4,i);hold on;
            mintime=min(Exp(i).time{1});
            for j=1:length(Exp(i).time)
                if Exp(i).JC==1
                    unidays=find(diff(Exp(i).time{j})>8); % JC time
                    unistart=[1 unidays+1];
                    uniend=[unidays length(Exp(i).time{j})];

                else
                    unistart=[];
                    uniend=[];
                    unidays=unique(floor(Exp(i).time{j})); % Ev time
                    for kk=1:length(unidays)
                        unistart(kk)=min(find(floor(Exp(i).time{j})==unidays(kk)));
                        uniend(kk)=max(find(floor(Exp(i).time{j})==unidays(kk)));
                    end
                end
                for ii=1:length(unistart)
                    times=Exp(i).time{j}(unistart(ii):uniend(ii));
                    rpts=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                    hamfilttimes=[];hamfiltrpts=[];
                    for jj=1:length(times)-lenham
                        hamfilttimes(jj)=mean(times(jj:jj+lenham-1).*h');
                        hamfiltrpts(jj)=mean(rpts(jj:jj+lenham-1).*h');
                    end
                        if Exp(i).JC==1
                            hamfilttimes=(hamfilttimes-mintime)/24;
                        else
                            hamfilttimes=hamfilttimes-mintime;
                        end
                    
                    if Exp(i).pre==j
                        plot(hamfilttimes,hamfiltrpts,'b.','Markersize',15)
                    else if Exp(i).wn==j
                            plot(hamfilttimes,hamfiltrpts,'r.','Markersize',15)
                        else if Exp(i).post==j
                                plot(hamfilttimes,hamfiltrpts,'b.','Markersize',15)
                            end
                        end
                    end
                end

            end
        end
     
        
        subplot(322);hold on;
        stairs(b,a,'Linewidth',3)
        stairs(d,c,'color','r','linewidth',3)
        %stairs(f,e,'color','k','linewidth',3)
        xlim([0 20])
        subplot(324);hold on;
        stairs(b,c-a,'color','r','Linewidth',3)
        plot([0 20],[0 0],'k')
        winwidth=50;
        subplot(326);plot(runningaverage(timevalsPRE1017-min(timevalsPRE1017),winwidth),runningaverage(distntPRE1017,winwidth),'*') % 1015_evtafreptest
        hold on;plot(runningaverage(timevalsWN1018(1:189)-min(timevalsPRE1017),winwidth),runningaverage(distntWN1018(1:189),winwidth),'*','Color','r') % 1017_wnonrepeats
        hold on;plot(runningaverage(timevalsWN1018(190:end)-min(timevalsPRE1017),winwidth),runningaverage(distntWN1018(190:end),winwidth),'*','Color','r') % 1017_wnonrepeats
        hold on;plot(runningaverage(timevalsPost1020-min(timevalsPRE1017),winwidth),runningaverage(distntPost1020,winwidth),'*','Color','k') % 1017_wnonrepeats   
        plot([0 100],[mean(distntPRE1017) mean(distntPRE1017)])
%%%%%%%%%%%%%%        
figure;hold on;
for i=7:8
for j=1:3
plot(runningaverage(Exp1(i).time{j}(),20),runningaverage(Exp1(i).rplength{j}(),20))
end
end        
%%%%%%%%
        Exp=Exp1;
        figure;hold on;
        lenham=11;
        h=hamming(lenham);
        h=h*(lenham/sum(h));

        for i=[9 11]%1:length(Exp) %- 2 is bk1bk29 - better example
            for j=1:length(Exp(i).time)
                if Exp(i).JC==1
                    unidays=find(diff(Exp(i).time{j})>8); % JC time
                    unistart=[1 unidays+1];
                    uniend=[unidays length(Exp(i).time{j})];
                else
                    unistart=[];
                    uniend=[];
                    unidays=unique(floor(Exp(i).time{j})); % Ev time
                    for kk=1:length(unidays)
                        unistart(kk)=min(find(floor(Exp(i).time{j})==unidays(kk)));
                        uniend(kk)=max(find(floor(Exp(i).time{j})==unidays(kk)));
                    end
                end
                
                    times=Exp(i).time{j}(unistart(ii):uniend(ii));
                    rpts=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                    hamfilttimes=[];hamfiltrpts=[];
                    for jj=1:length(times)-lenham
                        hamfilttimes(jj)=mean(times(jj:jj+lenham-1).*h');
                        hamfiltrpts(jj)=mean(rpts(jj:jj+lenham-1).*h');
                    end
                    if Exp(i).pre==j
                        plot(hamfilttimes-mintime,hamfiltrpts,'b.','Markersize',15)
                    else if Exp(i).wn==j
                            plot(hamfilttimes-mintime,hamfiltrpts,'r.','Markersize',15)
                        else if Exp(i).post==j
                                plot(hamfilttimes-mintime,hamfiltrpts,'b.','Markersize',15)
                            end
                        end
                    end
                end

            end
        end
     
%%%%%%%%%5
        
        
        
        
        
        
%%%%
% bk1bk29 - lesioned - October 2010 - DECREASE
        clear all
        load /cardinal3/SyntaxBirds/bk1bk29/Repeats1017.mat
        subplot(321);hold on;
        stairs(b,a,'Linewidth',3)
        stairs(d,c,'color','r','linewidth',3)
        %stairs(f,e,'color','k','linewidth',3)
        xlim([0 20])
        subplot(323);hold on;
        stairs(b,c-a,'color','r','Linewidth',3)
        plot([0 20],[0 0],'k')
        winwidth=50;
        subplot(325);plot(runningaverage(timevalsPRE1016-min(timevalsPRE1016),winwidth),runningaverage(distntPRE1016,winwidth),'*') % 1015_evtafreptest
        hold on;plot(runningaverage(timevalsWN1017(1:206)-min(timevalsPRE1016),winwidth),runningaverage(distntWN1017(1:206),winwidth),'*','Color','r') % 1017_wnonrepeats
        hold on;plot(runningaverage(timevalsWN1017(207:end)-min(timevalsPRE1016),winwidth),runningaverage(distntWN1017(207:end),winwidth),'*','Color','r') % 1017_wnonrepeats
        hold on;plot(runningaverage(timevalsPost1019-min(timevalsPRE1016),winwidth),runningaverage(distntPost1019,winwidth),'*','Color','k') % 1017_wnonrepeats
        plot([0 100],[mean(distntPRE1016) mean(distntPRE1016)])

% 1. bk67bk44 - lesioned - October 2010 - DECREASE   
% 2. bk1bk29 - lesioned - October 2010 - DECREASE
% 3. bk67bk44 - lesioned - October 2010 - INCREASE (same note!)
load /cardinal3/SyntaxBirds/bk67bk44/RepeatExpansion1023.mat 

% 4. w22pk23 - decrease but held for 20 days so recovery is slow
% 5. w45pk27 - increase - also held for a long time
% 6. w45pk27 - decrease (same note!) - recovery data is missing
% 7. w93pk62 - decrease - all data present
% 8. w93pk62 - increase (same note as 7) --- good example
% 9. w93pk62 - decrease - different note