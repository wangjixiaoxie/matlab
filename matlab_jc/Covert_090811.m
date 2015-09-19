% load /cardinal5/Covert040810.mat
% 121tmp.m


% EXAMPLES - bk76b
% 8.10.2011
clear all

load /bulbul2/CovertAnalysis/Covert070611.mat
            figure;hold on;
            subplot(131);hold on;
            w1=210:240;

            mnpre=mean(mean(ExpCTL.pitchPre(w1,:)));
            t1=mean(ExpCTL.timePre);
            p1=mean(mean(ExpCTL.pitchPre(w1,:)))-mnpre;
            s1=std(mean(ExpCTL.pitchPre(w1,:)))/sqrt(length(ExpCTL.pitchPre(w1,:)));
            t2=mean(ExpCTL.timeWN(397:end));
            p2=median(median(ExpCTL.pitchWN(w1,397:end)))-mnpre; % last hour
            s2=std(mean(ExpCTL.pitchWN(w1,397:end)))/sqrt(length(ExpCTL.pitchWN(w1,397:end)));
            t3=mean(ExpCTL.timePost(1:20));
            p3=mean(mean(ExpCTL.pitchPost(w1,1:20)))-mnpre; % first hour
            s3=std(mean(ExpCTL.pitchPost(w1,1:20)))/sqrt(length(ExpCTL.pitchPost(w1,1:20)));
            rat1=runningaverage(ExpCTL.timePre(1:end),50)-min(ExpCTL.timeWN);
            rat2=runningaverage(ExpCTL.timeWN(:),50)-min(ExpCTL.timeWN);
            rat3=runningaverage(ExpCTL.timePost(1:25),10)-min(ExpCTL.timeWN)-1;
            raf1=runningaverage(median(ExpCTL.pitchPre(w1,1:end))-mnpre,50);
            raf2=runningaverage(median(ExpCTL.pitchWN(w1,:))-mnpre,50);
            raf3=runningaverage(median(ExpCTL.pitchPost(w1,1:25))-mnpre,10);
            plot([rat1 rat2(1)],[raf1 raf2(1)])
            plot(rat2,raf2,'r')
            plot([rat2(end) rat3],[raf2(end) raf3])
%             plot(ExpCTL.timePre(1:1:end)-min(ExpCTL.timeWN),median(ExpCTL.pitchPre(w1,1:1:end))-mnpre,'.')   
%             plot(ExpCTL.timeWN(1:1:end)-min(ExpCTL.timeWN),median(ExpCTL.pitchWN(w1,1:1:end))-mnpre,'r.')
%             plot(ExpCTL.timePost(1:1:60)-min(ExpCTL.timeWN)-1,median(ExpCTL.pitchPost(w1,1:1:60))-mnpre,'.')
            %plot(runningaverage(ExpCTL.timeWN(1:5:end)-min(ExpCTL.timeWN),50),runningprctile(median(ExpCTL.pitchWN(w1,1:5:end))-mnpre,60,50))
            %plot(runningaverage(ExpCTL.timeWN(1:5:end)-min(ExpCTL.timeWN),50),runningprctile(median(ExpCTL.pitchWN(w1,1:5:end))-mnpre,40,50))
            %plot(runningaverage(ExpCTL.timePost-min(ExpCTL.timeWN),20),runningprctile(median(ExpCTL.pitchPost(w1,:))-mnpre,60,20))
            %plot(runningaverage(ExpCTL.timePost-min(ExpCTL.timeWN),20),runningprctile(median(ExpCTL.pitchPost(w1,:))-mnpre,40,20))

            %errorbar(x,y,e)
            plot(x,y,'Linewidth',3)
            ylim([-150 180])
            xlim([-5 8])
            plot([-5 8],[0 0])
            % figure;subplot(132);hold on;
            % mnpre=mean(mean(pitchNTpre(w1,:)));
            % w1=250:350;
            % t1=mean(timeNTpre);
            % p1=mean(mean(pitchNTpre(w1,:)))-mnpre;
            % s1=std(mean(pitchNTpre(w1,:)))/sqrt(length(pitchNTpre(w1,:)));
            % t2=mean(timeNTwn(280:end));
            % p2=median(median(pitchNTwn(w1,280:end)))-mnpre; % last hour
            % s2=std(median(pitchNTwn(w1,280:end)))/sqrt(length(pitchNTwn(w1,280:end)));
            % t3=mean(timeNTpost);
            % p3=median(median(pitchNTpost(w1,:)))-mnpre; % first hour
            % s3=std(mean(pitchNTpost(w1,:)))/sqrt(length(pitchNTpost(w1,:)));
            % x=[t1 t2 t3];y=[p1 p2 p3];e=[s1 s2 s3];
            % errorbar(x,y,e)
            % ylim([-20 60])
            subplot(132);hold on;
            w1=242:306;
            mnpre=mean(mean(ExpAPV.pitchACpre(w1,:)));
            mnapv=mean(mean(ExpAPV.pitchAPV(w1,:)));
            t1=mean(ExpAPV.timeAPV);
            p1=mean(mean(ExpAPV.pitchAPV(w1,:)))-mnapv;
            s1=std(mean(ExpAPV.pitchAPV(w1,:)))/sqrt(length(ExpAPV.pitchAPV(w1,:)));
            t2=mean(ExpAPV.timeAPVwn(end-20:end));
            p2=mean(mean(ExpAPV.pitchAPVwn(w1,34:end)))-mnapv; % last hour
            s2=std(mean(ExpAPV.pitchAPVwn(w1,34:end)))/sqrt(length(ExpAPV.pitchAPVwn(w1,34:end)));
            t3=mean(ExpAPV.timeACpost(1:160));
            p3=mean(mean(ExpAPV.pitchACpost(w1,1:160)))-mnpre; % first hour
            s3=std(mean(ExpAPV.pitchACpost(w1,1:160)))/sqrt(length(ExpAPV.pitchACpost(w1,1:160)));

            plot([runningaverage(ExpAPV.timeACpre(:),15)-min(ExpAPV.timeAPVwn)+3 runningaverage(ExpAPV.timeAPVwn(:),15)-min(ExpAPV.timeAPVwn) runningaverage(ExpAPV.timeACpost(1:160),50)-min(ExpAPV.timeAPVwn)-12],...
                [runningaverage(median(ExpAPV.pitchACpre(w1,:))-mnpre,15) runningaverage(median(ExpAPV.pitchAPVwn(w1,:))-mnapv,15) runningaverage(median(ExpAPV.pitchACpost(w1,1:160))-mnpre,50)])
             plot(ExpAPV.timeACpre(:)-min(ExpAPV.timeAPVwn)+3,median(ExpAPV.pitchACpre(w1,:))-mnpre,'.')
            plot(ExpAPV.timeAPVwn(:)-min(ExpAPV.timeAPVwn),median(ExpAPV.pitchAPVwn(w1,:))-mnapv,'r.')
            plot(ExpAPV.timeACpost(1:2:160)-min(ExpAPV.timeAPVwn)-12,median(ExpAPV.pitchACpost(w1,1:2:160))-mnpre,'.')

            %errorbar(x,y,e,'r')
            plot(x,y,'Linewidth',3)
            ylim([-100 150])
            xlim([-5 8])
            plot([-5 8],[0 0])
figure;
            subplot(133);hold on;
            w1=200:250;
            mnpre=mean(mean(ExpAPV.pitchACpreCTL(w1,:)));
            mnapv=mean(mean(ExpAPV.pitchAPVCTL(w1,:)));
            t1=mean(ExpAPV.timeAPV);
            p1=mean(mean(ExpAPV.pitchAPVCTL(w1,:)))-mnapv;
            s1=std(mean(ExpAPV.pitchAPVCTL(w1,:)))/sqrt(length(ExpAPV.pitchAPVCTL(w1,:)));
            t2=mean(ExpAPV.timeAPVwnCTL(34:end));
            p2=mean(mean(ExpAPV.pitchAPVwnCTL(w1,34:end)))-mnapv; % last hour
            s2=std(mean(ExpAPV.pitchAPVwnCTL(w1,34:end)))/sqrt(length(ExpAPV.pitchAPVwnCTL(w1,34:end)));
            t3=mean(ExpAPV.timeACpostCTL(1:30));
            p3=mean(mean(ExpAPV.pitchACpostCTL(w1,1:30)))-mnpre; % first hour
            s3=std(mean(ExpAPV.pitchACpostCTL(w1,1:30)))/sqrt(length(ExpAPV.pitchACpostCTL(w1,1:30)));
            rat1=runningaverage(ExpAPV.timeACpreCTL(1:end),40)-min(ExpAPV.timeAPVwnCTL)+2;
            rat2=runningaverage(ExpAPV.timeAPVwnCTL(:),20)-min(ExpAPV.timeAPVwnCTL);
            rat3=runningaverage(ExpAPV.timeACpostCTL(1:34),15)-min(ExpAPV.timeAPVwnCTL)-12;
            raf1=runningmedian(median(ExpAPV.pitchACpreCTL(w1,1:end))-mnpre,40);
            raf2=runningaverage(median(ExpAPV.pitchAPVwnCTL(w1,:))-mnapv,20);
            raf3=runningaverage(median(ExpAPV.pitchACpostCTL(w1,1:34))-mnpre,15);
            plot([rat1 rat2(1)],[raf1 raf2(1)])
            plot(rat2,raf2,'r')
            plot([rat2(end) rat3],[raf2(end) raf3])

%             plot([t1-t1 runningaverage(ExpAPV.timeAPVwnCTL(:),20)-min(ExpAPV.timeAPVwnCTL) runningaverage(ExpAPV.timeACpostCTL(1:34),15)-min(ExpAPV.timeAPVwnCTL)-12],[p1 runningaverage(median(ExpAPV.pitchAPVwnCTL(w1,:))-mnapv,20) runningaverage(median(ExpAPV.pitchACpostCTL(w1,1:34))-mnpre,15)])
%             plot(ExpAPV.timeAPVwnCTL(:)-min(ExpAPV.timeAPVwnCTL),median(ExpAPV.pitchAPVwnCTL(w1,:))-mnapv,'r.')
%             plot(ExpAPV.timeACpostCTL(1:34)-min(ExpAPV.timeAPVwnCTL)-12,median(ExpAPV.pitchACpostCTL(w1,1:34))-mnpre,'.')

            x=[t1-t1 t2-t1 t2+2-t1];y=[p1 p2 p3];e=[s1 s2 s3];
           % errorbar(x,y,e,'k')
           x=[0 4 6];
            plot(x,y,'Linewidth',3)
            ylim([-50 100])
            xlim([-5 8])
            plot([-5 8],[0 0])
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NONE OF THE ONES BELOW ARE VERY GOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
clear all
load /bulbul2/CovertAnalysis/Experiment14.mat
            figure;subplot(132);hold on;
            w1=290:355;
            mnpre=mean(mean(ExpAPV.pitchACpre(w1,:)));
            mnapv=mean(mean(ExpAPV.pitchAPV(w1,:)));
            t1=mean(ExpAPV.timeAPV);
            p1=mean(mean(ExpAPV.pitchAPV(w1,:)))-mnapv;
            s1=std(mean(ExpAPV.pitchAPV(w1,:)))/sqrt(length(ExpAPV.pitchAPV(w1,:)));
            t2=mean(ExpAPV.timeAPVwn(end-20:end));
            p2=mean(mean(ExpAPV.pitchAPVwn(w1,50:end)))-mnapv; % last hour
            s2=std(mean(ExpAPV.pitchAPVwn(w1,50:end)))/sqrt(length(ExpAPV.pitchAPVwn(w1,50:end)));
            t3=mean(ExpAPV.timeACpost(32:300));
            p3=mean(mean(ExpAPV.pitchACpost(w1,32:300)))-mnpre; % first hour
            s3=std(mean(ExpAPV.pitchACpost(w1,32:300)))/sqrt(length(ExpAPV.pitchACpost(w1,32:300)));
            plot([t1-t1 runningaverage(ExpAPV.timeAPVwn(:),40)-min(ExpAPV.timeAPVwn) t2-t1],[p1 runningaverage(median(ExpAPV.pitchAPVwn(w1,:))-mnapv,40) p2])
            plot([t2-t1 runningaverage(ExpAPV.timeACpost(32:300),80)-min(ExpAPV.timeAPVwn)],[p2 runningaverage(median(ExpAPV.pitchACpost(w1,32:300))-mnpre,80)])
            plot(ExpAPV.timeAPVwn(:)-min(ExpAPV.timeAPVwn),median(ExpAPV.pitchAPVwn(w1,:))-mnapv,'.')
            plot(ExpAPV.timeACpost(32:2:300)-min(ExpAPV.timeAPVwn),median(ExpAPV.pitchACpost(w1,32:2:300))-mnpre,'.')

            x=[t1-t1 t2-t1 t2+2-t1];y=[p1 p2 p3];e=[s1 s2 s3];
            x=[0 4 6];
            %errorbar(x,y,e,'r')
          figure;  plot(x,y,'Linewidth',3)
            ylim([-100 150])
            xlim([-2 10])
            plot([-2 10],[0 0])

clear all
load /bulbul2/CovertAnalysis/Exp7.mat %/cardinal6/CovertAnalysis/Exp7.mat
% NOT GOOD
ExpAPV=Experiment7;
            figure;subplot(132);hold on;
            w1=140:170;
            mnpre=mean(mean(ExpAPV.pitchACpre(w1,114:end)));
            mnapv=mean(mean(ExpAPV.pitchAPV(w1,:)));
            t1=mean(ExpAPV.timeAPV);
            p1=mean(mean(ExpAPV.pitchAPV(w1,:)))-mnapv;
            s1=std(mean(ExpAPV.pitchAPV(w1,:)))/sqrt(length(ExpAPV.pitchAPV(w1,:)));
            t2=mean(ExpAPV.timeAPVwn(end-20:end));
            p2=mean(mean(ExpAPV.pitchAPVwn(w1,75:end)))-mnapv; % last hour
            s2=std(mean(ExpAPV.pitchAPVwn(w1,75:end)))/sqrt(length(ExpAPV.pitchAPVwn(w1,75:end)));
            t3=mean(ExpAPV.timeACpost(1:160));
            p3=mean(mean(ExpAPV.pitchACpost(w1,1:160)))-mnpre; % first hour
            s3=std(mean(ExpAPV.pitchACpost(w1,1:160)))/sqrt(length(ExpAPV.pitchACpost(w1,1:160)));
            plot([t1-t1 runningaverage(ExpAPV.timeAPVwn(:),40)-min(ExpAPV.timeAPVwn) t2-t1],[p1 runningaverage(median(ExpAPV.pitchAPVwn(w1,:))-mnapv,40) p2])
            plot([t2-t1 runningaverage(ExpAPV.timeACpost(1:160),80)-min(ExpAPV.timeAPVwn)],[p2 runningaverage(median(ExpAPV.pitchACpost(w1,1:160))-mnpre,80)])
            plot(ExpAPV.timeAPVwn(:)-min(ExpAPV.timeAPVwn),median(ExpAPV.pitchAPVwn(w1,:))-mnapv,'.')
            plot(ExpAPV.timeACpost(1:2:160)-min(ExpAPV.timeAPVwn),median(ExpAPV.pitchACpost(w1,1:2:160))-mnpre,'.')

            x=[t1-t1 t2-t1 t2+2-t1];y=[p1 p2 p3];e=[s1 s2 s3];
            x=[0 4 6];
            %errorbar(x,y,e,'r')
            plot(x,y,'Linewidth',3)
            ylim([-100 150])
            xlim([-2 10])
            plot([-2 10],[0 0])

figure;hold on;
ravgwin1=40;    window1=[140:170];
% Targeted syllable - mean FF            % THR ~= 2200Hz
subplot(411);hold on;   plotExp7targmean(Experiment7,ravgwin1,window1); ylim([2000 2400])
% Targeted syllable - s.d. of FF        (after removing outliers)
subplot(412);hold on;   plotExp7targsd(Experiment7,ravgwin1,window1); ylim([35 90])
% Non-targeted syllable - mean FF
window1=[170:340];
subplot(413);hold on;   plotExp7nontargmean(Experiment7,ravgwin1,window1); ylim([1900 2300])



% Example experiment - Positive control learning w/o APV [bk76bk63 - ExperimentPC(2)].
% This example is nice for showing rapid smooth FF learning.  Only downside is
% rapid decay of learning after WN off.
        load /cardinal6/bk76bk63/EPC2.mat
        figure;hold on;
            window1=160:240;
            ravgwin=20;
            plotExpPC2(EPC2,ravgwin,window1);plot([8004.4 8008.3],[2510 2510])
            xlim([8002 8012.65]);ylim([2350 2650])
            
% plot all PCs - also see ExpPC2 and b4o59
        figure;hold on;
        for i=1:14
            subplot(4,4,i);hold on;
            if length(ExperimentPC(i).time)==1
            window1=ExperimentPC(i).time-30:ExperimentPC(i).time;
            else
                window1=ExperimentPC(i).time;
            end
            plot(runningmedian(ExperimentPC(i).timePre,10),runningmedian(mean(ExperimentPC(i).pitchPre(window1,:)),10),'b.')
            plot(runningmedian(ExperimentPC(i).timeWN,10),runningmedian(mean(ExperimentPC(i).pitchWN(window1,:)),10),'r.')
              plot(runningmedian(ExperimentPC(i).timePost,10),runningmedian(mean(ExperimentPC(i).pitchPost(window1,:)),10),'k.')
        end
            
% Example experiment (b10o7 - #7) showing how any acute effects on FF go
% away quickly in non-targeted syllable, but learning persists in targeted
% sylalble.  Also shows reduction of FF variation (in targeted syllable).
        load /bulbul2/CovertAnalysis/Exp7.mat %/cardinal6/CovertAnalysis/Exp7.mat
        figure;hold on;
        ravgwin1=40;    window1=[140:170];
        % Targeted syllable - mean FF            % THR ~= 2200Hz
            subplot(411);hold on;   plotExp7targmean(Experiment7,ravgwin1,window1); ylim([2000 2400])           
        % Targeted syllable - s.d. of FF        (after removing outliers)   
            subplot(412);hold on;   plotExp7targsd(Experiment7,ravgwin1,window1); ylim([35 90])
        % Non-targeted syllable - mean FF           
            window1=[170:340];
            subplot(413);hold on;   plotExp7nontargmean(Experiment7,ravgwin1,window1); ylim([1900 2300])
        % Non-targeted syllable - s.d. of FF - crappy  
            subplot(414);hold on;   plotExp7nontargsd(Experiment7,ravgwin1,window1); ylim([35 110])
            
            
            
            
            
            
            
            
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
load /bulbul4/CovertAnalysis/CovertCompact4.mat
% bar graphs
    figure;hold on;
    bar(0.5,100*mean(finalpointPC))
    plot(0.5,100*(finalpointPC),'.','markersize',15)
    bar(1.5,100*mean(finalpoint(ind)))
    plot(1.5,100*(finalpoint(ind)),'.','markersize',15)  
    bar(3,100*mean(postPC1))
    plot(3,100*(postPC1),'.','markersize',15)  
    bar(4,100*mean(post2(ind)))
    plot(4,100*(post2(ind)),'.','markersize',15)  
    errorbar(0.5,100*mean(finalpointPC),100*semPCfinal)
    errorbar(1.5,100*mean(finalpoint(ind)),100*semfinal)
    errorbar(3,100*mean(postPC1),100*(std(postPC1)/sqrt(length(postPC1))))
    errorbar(4,100*mean(post2(ind)),100*sempost)
    xlim([0 4.5])
    ylim([-2 4])
% 
    figure;hold on;
    bar(0.5,100*mean(post2(ind)))
    plot(0.5,100*(post2(ind)),'.','markersize',15)  
    bar(1.5,100*mean(Post2NT))
    plot(1.5,100*(Post2NT),'.','markersize',15)  
    errorbar(0.5,100*mean(post2(ind)),100*sempost)
    errorbar(1.5,100*mean(Post2NT),100*semNT)    
    xlim([0 2])
    ylim([-0.5 1.5])


% Krupa figure with confidence intervals
% In this representation, error bars are resampled CIs, time point is
    % entire second day, positive control time point is last 15 syllables
    % during wn.
figure;hold on;
                % not using CIs
plotbasicKrupNTnCI2(Post2NT,ind,finalpointPC,finalpoint,post2,semPCfinal,semfinal,sempost,semNT,pcCIbase,CIbaseCTL,CIbase,postPC1)
xlim([-0.5 4])
    %plotbasicKrup(ind,pcCIbase,pcCIfinal,finalpointPC,CIbase,CIfinal,CIpost1,finalpoint,post1,CIpost2,pcCIpost1,pcCIpost2)  
   figure;hold on; plotbasicKrupNT(CIbaseCTL,CIpost1NT,CIpost2NT,Post1NT,Post2NT,ind,pcCIbase,pcCIfinal,finalpointPC,CIbase,CIfinal,CIpost1,finalpoint,post1,CIpost2,pcCIpost1,pcCIpost2)  
% In this representation, resampled CIs, time point 1 is last 50 of first
    % day (for exps when there is data), time point 2 is first 50 of second day
    % (for all exps).
    
    figure;hold on;
    %plotsleepKrup(CIpostday1,CIpostday2,mnpostday1,mnpostday2,ind,pcCIbase,pcCIfinal,finalpointPC,CIbase,CIfinal,CIpost1,finalpoint,post1,CIpost2,pcCIpost1,pcCIpost2)  
    plotsleepKrupNT(CIpost1NT,CIpost2NT,Post1NT,Post2NT,CIpostday1,CIpostday2,mnpostday1,mnpostday2,ind,pcCIbase,pcCIfinal,finalpointPC,CIbase,CIfinal,CIpost1,finalpoint,post1,CIpost2,pcCIpost1,pcCIpost2)      
% bar graphs 
    % Block of expression: finalpointPC vs. finalpoint, ttest2 p-value is 0.0100, ranksum p-value is 0.0148
    % Block of acquisition: finalpointPC vs. mnpostday2, ttest2 p-value is 0.8104, ranksum p-value is 0.9597
    % finalpoint vs. mnpostday2, ttest2 p-value is 0.0073, ranksum p-value is 0.0037
    figure;hold on;bar([1 2 3 4],[mean(finalpointPC) mean(finalpoint(ind)) mean(post2(ind)) mean(Post2NT)])
    plot([1 1],[mean(finalpointPC)-semPCfinal mean(finalpointPC)+semPCfinal],'r');plot([2 2],[mean(finalpoint(ind))-semfinal mean(finalpoint(ind))+semfinal],'r');plot([3 3],[mean(post2(ind))-sempost mean(post2(ind))+sempost],'r');plot([4 4],[mean(Post2NT)-semNT mean(Post2NT)+semNT],'r')
    plot(1,finalpointPC,'b.','Markersize',15);plot(2,finalpoint(ind),'b.','Markersize',15);plot(3,post2(ind),'b.','Markersize',15);plot(4,Post2NT,'b.','Markersize',15)
            % postday - index of beginning and end of next day for targeted notes
            % postday2 - index of beginning and end of next day for control notes
            % ind - all experiments
            % ind2 - experiments where we have data for control notes
            % ind3
            % LearnCoef - equal to +1 if adaptive direction is upwards
            % FFmwin - measurement window for targeted notes
            % FFmwinCTL - measurement window for CTL notes
            % FFacpreCTL - mean FF in the measurement window
            % FFacpostCTL - mean FF in the measurement window
            % FFapvCTL - mean FF in the measurement window'
            % mpitchpreTarg - mean baseline pitch of targeted note - handy for converted from Hz to % and vice versa
            % mpitchprePC - mean baseline pitch of positive control notes
            % mpitchpreNT - mean baseline pitch of non-targeted notes
            % post1 - LEARNING at some early time point (first 50 after apv off, I think)
            % post2 - LEARNING measured on the entire 2nd day
            
% Variability reduction in targeted notes
            % Variability in CTL notes
                % VarPreCTL - CV in control notes
                % VarAPVCTL - CV in control notes (with APV)
                % VarPost1CTL - CV of final 50 notes of the day when APV+wn ended - if there are less than 50 (e.g. zero) then this variable is empty
                % VarPost2CTL - CV of first 50 notes of the day after APV+wn ended
                % DayoneCTL - includes raw data that goes into VarPost1CTL
                % DaytwoCTL - includes raw data that goes into VarPost2CTL
            % Variability in targeted notes
                % VarPre - CV 
                % VarAPV - CV  (with APV)
                % VarPost1 - CV of final 50 notes of the day when APV+wn ended - if there are less than 50 (e.g. zero) then this variable is empty
                % VarPost2 - CV of first 50 notes of the day after APV+wn ended
                % Dayone - includes raw data that goes into VarPost1
                % Daytwo - includes raw data that goes into VarPost2
plotVarPoints(Experiment7,VarPre,VarPost,VarAPV,VarPost1,VarPost2,VarPreCTL,VarPostCTL,VarAPVCTL,VarPost1CTL,VarPost2CTL,ind,ind2)
plotVarLines(...) % need to include sem bars
        
    % CONCLUSION - variability pretty much recovers by the end of the first day
        
            mean(VarAPV(ind)./VarPre(ind))
            % Var reduction in control (non-targeted notes)?
            mean(VarAPVCTL(ind2)./VarPreCTL(ind2));bbb=VarPost2CTL(ind2)./VarPreCTL(ind2);
            mean(bbb(~isnan(bbb))) % 0.9645

% Acute vs. lasting FF changes in non-targeted control notes
            % Acute - acute change in pitch of CTL notes (diff b/t APV and pre)
            % FFpermanent1 - shift in CTL notes measured at end of day (diff b/t post1 and pre)
            % Permanent - lasting shift in CTL notes measured on 2nd day (diff b/t post and pre)
            % LearningCTL - Permanent*LearnCoef (is permanent shift in the adaptive direction) - addresses specificity of learning
            % AcuteTarg - acute change in pitch of targeted notes
            % FFafter1targ - in targeted notes, mean FF of last 50 notes of the day when APV+wn ended - if there are notes in that time window
            % FFafter2targ - in targeted notes, mean FF of first 50 notes of the day after APV+wn ended
            % FFbeforetarg - in targeted notes, mean FF pre
        % No relationship between Acute APV-dep FF offset in CTL note and Permanent FF offset in CTL note
            corrcoef(Acute(ind2),Permanent(ind2)) % R = 0.0852
            corrcoef(Acute(find(FFpermanent1~=0)),FFpermanent1(find(FFpermanent1~=0))) % R = 0.16 (n=7) - no relationship   
figure;hold on;
plot(1,Acute(ind2)./mpitchpreNT(ind2),'b.','Markersize',15);plot(2,Permanent(ind2)./mpitchpreNT(ind2),'b.','Markersize',15)  
plot(1,mean(Acute(ind2)./mpitchpreNT(ind2)),'b.','Markersize',35);plot(2,mean(Permanent(ind2)./mpitchpreNT(ind2)),'b.','Markersize',35)
plot(1,mean(abs(Acute(ind2))./mpitchpreNT(ind2)),'r.','Markersize',35);plot(2,mean(abs(Permanent(ind2))./mpitchpreNT(ind2)),'r.','Markersize',35)
plot([1;2],[((Acute(ind2))./mpitchpreNT(ind2));((Permanent(ind2))./mpitchpreNT(ind2))],'b')
plot([0 3],[0 0],'k')
xlim([0.5 2.5])
%%%%%%%%%%%%
%%%%%%%%%%%%%
figure;hold on;
plotPreSleepFigure(Permanent,VarAPV,VarAPVCTL,indDN,Learn1,Learn2,mpitchpreTarg,ind3,Acute,FFpermanent1,VarPreCTL,VarPost1CTL,VarPost1,VarPost2CTL,VarPost2,VarPre)
% CONCLUSION - acute effect totally recovers by end of first day
                        gq=(FFpermanent1./mpitchpreNT).*LearnCoef;
                        Post1NT=gq([1 2 7 10 13 15 16]);
                        gq=(Permanent./mpitchpreNT).*LearnCoef;
                        Post2NT=gq(ind2);
                        mean(Post1NT) % -0.0051 as proportion = -14.0Hz   % use this in krupa?
                        mean(Post2NT) % -0.0014 as proportion = -4.3 Hz   % use this in krupa?
            %          % 
            %             figure;hold on;plot(Acute(ind2),Permanent(ind2),'b.');plot([0 0],[-60 80]);plot([-40 80],[0 0]);plot(mean(Acute(ind2)),mean(Permanent(ind2)),'b+','Markersize',25)
            %             plot(Acute(find(FFpermanent1~=0)),FFpermanent1(find(FFpermanent1~=0)),'g.');plot([0 0],[-60 80]);plot([-40 80],[0 0]);plot(mean(Acute(find(FFpermanent1~=0))),mean(FFpermanent1(find(FFpermanent1~=0))),'g+','Markersize',25)
            %         % CTL notes don't shift in the learning direction (control for general shift of song FF)
            %             corrcoef(LearningCTL(ind2),Learning(ind2)) % R = 0.0038
            %             figure;plot(LearningCTL(ind2),Learning(ind2),'.','Markersize',15);hold on;plot([-80 40],[0 0]);plot([0 0],[-80 100])
            
            
% Some sort of decay representation
    figure;hold on;bar([1 2],[mean(LearnEnd(indNext)./mpitchprePC(indNext)) mean(mnNext(indNext)./mpitchprePC(indNext))])
    semPCday1=std(LearnEnd(indNext)./mpitchprePC(indNext))/sqrt(10);semPCday2=std(mnNext(indNext)./mpitchprePC(indNext))/sqrt(10);
    plot([1 1],[mean(LearnEnd(indNext)./mpitchprePC(indNext))-semPCday1 mean(LearnEnd(indNext)./mpitchprePC(indNext))+semPCday1],'r');plot([2 2],[mean(mnNext(indNext)./mpitchprePC(indNext))-semPCday2 mean(mnNext(indNext)./mpitchprePC(indNext))+semPCday2],'r');    
    plot(1,LearnEnd(indNext)./mpitchprePC(indNext),'b.','Markersize',15);plot(2,mnNext(indNext)./mpitchprePC(indNext),'b.','Markersize',15)
    plot([1;2],[LearnEnd(indNext)./mpitchprePC(indNext);mnNext(indNext)./mpitchprePC(indNext)],'b')
        % mnNext is the average FF for the entire "next day" after the end of
            % WN for the 10 experiments (indNext) for which I recorded song on the
            % next day
        % LearnEnd is the average FF for the last 15 songs with WN on
    figure;plot(LearnEnd(indNext)./mpitchprePC(indNext),mnNext(indNext)./mpitchprePC(indNext),'.')
    hold on;plot([-0.02 0.04],[-0.02 0.04],'k')
    
% Some sort of localization representation - really no localization
    for i=1:size(centered,1)
        indtarg=find(centered(i,490:510)~=0);
        mntarg(i)=mean(centered(i,indtarg+489));
        indall=find(centered(i,400:600)~=0);
        mncent(i)=mean(centered(i,indall+399));
    end
    for i=1:size(centeredPC,1)
        indtargPC=find(centeredPC(i,490:510)~=0);    
        mntargPC(i)=mean(centeredPC(i,indtargPC+489));    
        indallPC=find(centeredPC(i,400:600)~=0);
        mncentPC(i)=mean(centeredPC(i,indallPC+399));
    end
