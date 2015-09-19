% EXAMPLES WITH APV

            load /bulbul/cardinal8/Covert/CovertExample_060411.mat
            % Experiment 12 - bk75bk62 - 9.27.09
                % note that APVwn are catch trials - everything else is not
                % Threshold set at 

                
                
                
                
                
            figure;hold on;
            mrksize=20;
                subplot(211);hold on;
                coef=-1;
                m1=round(median(Exp12.TargetingWN))-64:round(median(Exp12.TargetingWN));
                mnbase=mean(mean(Exp12.pitchACpre(m1,:)));
                %m1=260:290;
                ravgwin=50;
                % All acsf trials from that morning
                plot(timing4(fvalsACpre),coef*(mean(pitchACpre(m1,:))-mnbase),'.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsACpre),ravgwin),runningaverage(coef*(mean(pitchACpre(m1,:))-mnbase),ravgwin),'k')
figure;hold on;
% plot
% plot
% plot
mnbase=mean(mean(pitchACpre(m1,100:end)));
% errorbar(mean(timing4(fvalsACpre)),0,std(coef*(mean(pitchACpre(m1,:))-mnbase))/sqrt(length(pitchACpre(m1,:))))       
% errorbar(mean(timing4(fvalsAPV)),mean(coef*(mean(pitchAPV(m1,:))-mnbase)),std(coef*(mean(pitchAPV(m1,:))-mnbase))/sqrt(length(pitchAPV(m1,:))))       
% errorbar(mean(timing4(fvalsAPVwnCatch(52:end))),mean(coef*(mean(pitchAPVwnCatch(m1,52:end))-mnbase)),std(coef*(mean(pitchAPVwnCatch(m1,52:end))-mnbase))/sqrt(length(pitchAPVwnCatch(m1,52:end))))   
% errorbar(mean(timing4(fvalsACpost(34:146))),mean(coef*(mean(pitchACpost(m1,34:146))-mnbase)),std(coef*(mean(pitchACpost(m1,34:146))-mnbase))/sqrt(length(pitchACpost(m1,34:146))))       

plot(mean(timing4(fvalsACpre)),0,'.','Markersize',25)       
plot(mean(timing4(fvalsAPV(1:100))),mean(coef*(mean(pitchAPV(m1,1:100))-mnbase)),'.','Markersize',25)       
plot(mean(timing4(fvalsAPVwnCatch(52:end))),mean(coef*(mean(pitchAPVwnCatch(m1,52:end))-mnbase)),'.','Markersize',25)       

plot(mean(timing4(fvalsACpost(34:146))),mean(coef*(mean(pitchACpost(m1,34:146))-mnbase)),'.','Markersize',25)       

                
                % Excluding first 30 minutes of apv - all notes sung
                plot(timing4(fvalsAPV),coef*(mean(pitchAPV(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsAPV),ravgwin),runningaverage(coef*(mean(pitchAPV(m1,:))-mnbase),ravgwin),'k')
                % All notes sung with apv + wn
                plot(timing4(fvalsAPVwnCatch),coef*(mean(pitchAPVwnCatch(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsAPVwnCatch),ravgwin),runningaverage(coef*(mean(pitchAPVwnCatch(m1,:))-mnbase),ravgwin),'k')
                % All acsf trials from the next morning
                plot(timing4(fvalsACpost(34:146)),coef*(mean(pitchACpost(m1,34:146))-mnbase),'.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsACpost(34:146)),ravgwin),runningaverage(coef*(mean(pitchACpost(m1,34:146))-mnbase),ravgwin),'k')
                plot([6505 6515],[coef*mean(mean(Exp12.pitchACpost(m1,34:146))-mnbase) coef*mean(mean(Exp12.pitchACpost(m1,34:146))-mnbase)],'k')    
                plot([6496 6500],[thrsh-mnbase thrsh-mnbase])
                plot([6485 6515],[0 0])
                ylim([-200 200])
                xlim([6485 6515])

                time_offset=min(timing4(fvalsAPV))-min(Exp12.timeAPVCTL);
                subplot(212);
                hold on;
                coef=-1;
                m1=250:300;
                pitchACpreCTL=pitchACpreCTL_1(:,[1:4:660]);
                timeACpreCTL=timing4(fvalsACpreCTL_1([1:4:660]));
                mnbase=mean(mean(pitchACpreCTL(m1,:)));
                %m1=round(median(Exp12.TargetingWN))-64:round(median(Exp12.TargetingWN));
                %m1=250:350;
                ravgwin=20;
                % All acsf trials from that morning
                plot(time_offset+timeACpreCTL+864.1,coef*(mean(pitchACpreCTL(m1,:))-mnbase),'.','Markersize',mrksize)
                plot(time_offset+runningaverage(timeACpreCTL+864.1,ravgwin),runningaverage(coef*(mean(pitchACpreCTL(m1,:))-mnbase),ravgwin),'k')
                % Excluding first 30 minutes of apv - all notes sung
                plot(time_offset+Exp12.timeAPVCTL,coef*(mean(Exp12.pitchAPVCTL(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(time_offset+runningaverage(Exp12.timeAPVCTL,ravgwin),runningaverage(coef*(mean(Exp12.pitchAPVCTL(m1,:))-mnbase),ravgwin),'k')
                % All notes sung with apv + wn
                plot(time_offset+Exp12.timeAPVwnCTL,coef*(mean(Exp12.pitchAPVwnCTL(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(time_offset+runningaverage(Exp12.timeAPVwnCTL,ravgwin),runningaverage(coef*(mean(Exp12.pitchAPVwnCTL(m1,:))-mnbase),ravgwin),'k')
                % All acsf trials from the next morning
                plot(time_offset+Exp12.timeACpostCTL(33:142),coef*(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase),'.','Markersize',mrksize)
                plot(time_offset+runningaverage(Exp12.timeACpostCTL(33:142),ravgwin),runningaverage(coef*(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase),ravgwin),'k')
                plot([6485 6515],[0 0])
                plot([6505 6515],[coef*mean(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase) coef*mean(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase)],'k')
                ylim([-200 200])
                xlim([6485 6515])

   % EXAMPLES W/O APV
   load /cardinal8/CovertNTExample_060411.mat
   close all
   for i=1:length(ExperimentPC)
       figure;hold on;
       coef=1-2*isequal(ExperimentPC(i).DIR,'down');
       if length(ExperimentPC(i).time)==1
           time=ExperimentPC(i).time;
       else
           time=round(mean(ExperimentPC(i).time));
       end
       mnbase=mean(mean(ExperimentPC(i).pitchPre(time,:)));
       plot(runningaverage(ExperimentPC(i).timePre,5),coef*runningmedian(ExperimentPC(i).pitchPre(time,:)-mnbase,5),'.')
       plot(runningaverage(ExperimentPC(i).timeWN,5),coef*runningmedian(ExperimentPC(i).pitchWN(time,:)-mnbase,5),'r.')       
       plot(runningaverage(ExperimentPC(i).timePost,5),coef*runningmedian(ExperimentPC(i).pitchPost(time,:)-mnbase,5),'.')

   end
   % numbers 2,3,6,7 most promising
   load /cardinal8/APVinRA/cardinal6birds/bk76bk63/EPC2_60511.mat
   % EXPERIMENT 2
   % 11.30.09 - targeted note is C - first flat stack - auto-labeled 'x'
           i=2;
           pitchWN=EPC2.pitchWN(:,[1:4:504]);
           timeWN=(EPC2.timeWN([1:4:504]));

           figure;hold on;
           coef=1-2*isequal(EPC2.DIR,'down');
           if length(EPC2.time)==1
               time=EPC2.time;
           else
               time=round(mean(EPC2.time));
           end
           mnbase=mean(mean(EPC2.pitchPre(time,:)));
           ravgwin=50;
           plot((EPC2.timePre),coef*(EPC2.pitchPre(time,:)-mnbase),'.')
           plot((timeWN),coef*(pitchWN(time,:)-mnbase),'r.')
           plot((EPC2.timePost(1:82)),coef*(EPC2.pitchPost(time,1:82)-mnbase),'.')

           plot(runningaverage(EPC2.timePre,ravgwin),coef*runningmedian(EPC2.pitchPre(time,:)-mnbase,ravgwin),'k','Linewidth',3)
           plot(runningaverage(timeWN,ravgwin),coef*runningmedian(pitchWN(time,:)-mnbase,ravgwin),'k','Linewidth',3)
           %plot(runningaverage(EPC2.timePost(1:82),ravgwin),coef*runningmedian(EPC2.pitchPost(time,1:82)-mnbase,ravgwin),'k','Linewidth',3)
           plot([8009 8012],[mean(EPC2.pitchPost(time,1:82)-mnbase) mean(EPC2.pitchPost(time,1:82)-mnbase)])
           plot([8002 8012],[0 0])
           ylim([-200 200])
           xlim([8002 8012])
    % CTL note (nt = non-targeted) is the note after 'a' - second high
    % stack
           figure;hold on;
           coef=1-2*isequal(EPC2.DIR,'down');
           time=[250:350];
           mnbase=mean(mean(pitchPREnt(time,:)));
           ravgwin=50;
           plot(timing4(fvalsPREnt),coef*(mean(pitchPREnt(time,:))-mnbase),'.')
           plot(timing4(fvalsWNnt),coef*(mean(pitchWNnt(time,:))-mnbase),'r.')
           plot(timing4(fvalsPOSTnt),coef*(mean(pitchPOSTnt(time,:)-mnbase),'.')

           plot(runningaverage(EPC2.timePre,ravgwin),coef*runningmedian(EPC2.pitchPre(time,:)-mnbase,ravgwin),'k','Linewidth',3)
           plot(runningaverage(timeWN,ravgwin),coef*runningmedian(pitchWN(time,:)-mnbase,ravgwin),'k','Linewidth',3)
           %plot(runningaverage(EPC2.timePost(1:82),ravgwin),coef*runningmedian(EPC2.pitchPost(time,1:82)-mnbase,ravgwin),'k','Linewidth',3)
           plot([8009 8012],[mean(EPC2.pitchPost(time,1:82)-mnbase) mean(EPC2.pitchPost(time,1:82)-mnbase)])
           plot([8002 8012],[0 0])
           ylim([-200 200])
           xlim([8002 8012])