% b4o59
% possible stim candidate
% very massive bird
cleandir4('batch',50000,500,5,5);
mk_rmdata('batch.dcrd',1)
!csh yes | ./rmdata
%
load /cardinal/b4o59/templates1217.mat
mk_tempf('batch.keep',templaWN,2,'obs0');
get_trigt2('batch.keep',cntrngWN,0.4,128,1,1);
label_trigs('batch.keep','b','obs0',5000,1,1,5,30);

% segment at 10000 (default)

% templaStim + 0.4ms refractory period works well
    % recognizes 150-200ms before targeted note
    % first template is BTAF off of top of downsweep
    % second template is wacky stuff at bottom of downsweep
    % even better stuff could be done with an amplitude threshold 

% templaWN + 0.4ms refractory period works well
    % synshifts
    % s.d. of 4ms
    % hits 100%

% launchpad/data/jcharles/b4o59

%%%% TEST FOR FF LEARNING RATE
    load Data1221.mat
    figure;hold on;
    plot(runningaverage(valsPRE(:,1),10),runningmedian(valsPRE(:,2),10),'b*')
    plot(runningaverage(valsWN(:,1),10),runningmedian(valsWN(:,2),10),'r*')
    plot(runningaverage(valsPOST(:,1),10),runningmedian(valsPOST(:,2),10),'k*')

    % 12.20.10
        % Hit above 2294Hz at noon (median=2309Hz-0.5*31Hz)
    % 12.21.10
        % WN off at 2pm
    % Summary
        % Didn't sing almost anything on 12.20 (day one)
        % Learned rather quickly on 12.21 (day two) once he started singing
        % Try again after Christmas to make sure he gets used to WN
 fvals0106pre=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);   
 pitch0106pre=pitchcontour(fvals0106pre,2000,3000);
 vals0106pre=evtaf_freq('batch.catch',[2000,3000],'b',128,'obs0',1,1);
 
    % 01.05.11 - back from winter break, start screening again
    % 01.06.11 - median is now 2305Hz 
    %          - 2:55pm - wn on 100% to make sure he'll sing through - he did
    %          - 6:00pm - wn off
    % 01.07.11 - morning median = 2325Hz
    %          - noon - hit above 2310Hz = (2325Hz-0.5*31Hz)
    %          - 7pm - hit above 2285Hz
    % 01.08.11
    %          - 2pm - hit above 2265Hz
    % 01.09.11
    %          - 11:30am - wn off
    
    % Covert positive control
    % 01.18.11 - acquistion starting around 1:30pm
    % plan to do 4hr control on 1.19
    % 01.19.11
        % 1:15pm to 5:20pm - WN on - hit above 2340Hz - the experiment went
        % very well - smooth learning, ~20Hz, rapid decay afterwards
    
    % 01.21.11
        % noon - WN on - hit above 2350Hz; 1pm hit above 2340Hz
        % 8pm - WN off/lights off
    % 01.23.11 - dawn - WN on - hit below 2350Hz
        % 2pm - WN off
    % 01.24.11 - dawn - WN on - hit below 2350Hz (same as 1.23)
        % 2pm - WN off (same as 1.23)
    % 01.25.11 - 4pm - WN on - hit below 2360Hz
    % 01.26.11 - 10am - WN off
    
    % 02.03.11 - surgery - bilateral stim arrays in LMAN
    % 02.10.11 - leads in - transferred to stim rig with single stimulator
    % 02.14.11 - singing
    % 02.16.11 - decent template (but should be improved)
    % 02.17.11
    %   10am - Stim on - R23,L23,40ms del, 80ms dur, 150uA single stimulator - nothing
    %   1pm - Stim    - R23,L23,60ms del,90ms dur, 250uA single stimulator - nothing
    %   7:20pm - Stim - R14,L14,50ms del, same params
    
    % 02.18.11
    %   11am - Stim - R14,L14,30msdel,90ms dur, 400uA single stimulator
    %       around 6Volts
    %       seems to work fairly nicely - upward offset and sd reduction
    
    % 05.06.11
    %   plugged in & singing
    % 05.07.11
    % Made stim template with excellent targeting
    %   2:40pm - Stim on - R23,L23,20ms del, 60ms dur, 50uA per stimulator
    %   No effect
    
    % 05.08.11
    %  2pm - Stim on - R14,L14,20ms del, 60ms dur, 100uA per stimulator
    % Very nice effect
    %  3pm - 
    
    % 05.13.11
    % Great targeting for stim and wn
    % 11:45pm - wn on, 100% stim - hit below 2350Hz
    % 1pm - hit below 2360Hz - 
    
    
    
   figure;hold on;
% plot(runningmedian(valsPre121(:,1),100),runningmedian(valsPre121(:,2),100),'b.')
% plot(runningmedian(valsWN121(:,1),100),runningmedian(valsWN121(:,2),100),'r.')
plot(runningmedian(valsPost122(:,1),100),runningmedian(valsPost122(:,2),100),'k.')
plot(runningmedian(valsWN123(:,1),100),runningmedian(valsWN123(:,2),100),'r.')
plot(runningmedian(valsPost123(:,1),100),runningmedian(valsPost123(:,2),100),'k.')
plot(runningmedian(valsWN124(:,1),100),runningmedian(valsWN124(:,2),100),'r.')   
plot(runningmedian(vals125pre(:,1),100),runningmedian(vals125pre(:,2),100),'k.')
plot(runningmedian(vals126WN(:,1),100),runningmedian(vals126WN(:,2),100),'r.')   
plot(runningmedian(vals126post(:,1),100),runningmedian(vals126post(:,2),100),'k.')

    figure;plot(runningmedian(valsWN(:,1),50),runningmedian(valsWN(:,2),50),'r.')
hold on;plot(runningmedian(valsWN(:,1),50),runningmedian(valsWN(:,2),50),'r.')
hold on;plot(runningmedian(valsPost(:,1),50),runningmedian(valsPost(:,2),50),'k.')

% RECORDING

% CH0 - obs0 - song
% CH1 - obs1 -  R3 - open band (1Hz to 10kHz)
% CH2 - obs2 -  L3 - filtered (300Hz to 10kHz)
% CH3 - obs3 -  L4 - open band (1Hz to 10kHz)
% CH4 - obs4 -  R1 - filtered (300Hz to 10kHz)
% REF - silver wire in bird brain, buffered

% 04.17.11 - leads on
fvalsBOS=findwnoteSPK('batchnotes','x','','',0,[2000 2700],3e5,1,'obs0',0,1500);
BOSstim=fvalsBOS(1).datt(1:290000);
wavwrite(BOSstim/(max(abs(BOSstim))*1.0001),32000,16,'b4o59_BOS.wav')
% 04.20.11
    % 10am - put female bird in the cage
    % 12:30pm - put male bird in the cage
    % 2:30pm - solo
    % STILL NOT SINGING
    % 3:24pm (begin neural recording) - solo - playback WN, playback r30 song
    % 4:40pm - neural recording - playback BOS
% 04.21.11 - SINGING
    % made template for 'a' - not great, because3
    % 2:30pm - wn on 50% of renditions of syllable 'a'
% 04.22.11
    % 10:20am - turned off wn
    % 10:20am - CHANGED REF to L1
% CH0 - obs0 - song
% CH1 - obs1 -  R3 - open band (1Hz to 10kHz)
% CH2 - obs2 -  L3 - filtered (300Hz to 10kHz)
% CH3 - obs3 -  L4 - open band (1Hz to 10kHz)
% CH4 - obs4 -  R1 - filtered (300Hz to 10kHz)
% REF - L1

% 04.23.11
% 11am - tried R2 as ref briefly
% 2pm, REF = L1, hit below 2360Hz (70th percentile)
    

% 04.26.11
pitch426=pitchcontour(fv,6600,7500);
ptvls426=median(pitch426(880:1000,:))
figure;plot(ptvls426)
ind1=find(ptvls426>6800);
for i=1:size(Datt426(2).data,1)
gg2(i)=prctile(abs(Datt426(2).data(i,27000:33000)),90)-median(abs(Datt426(2).data(i,1:27000)));
end
hold on;plot(50*(runningaverage(gg2,50)-mean(gg2)),'r')

%%%%
%%%%


% 1. label song 'a'
% 2. run this code to pull out data from all channels
    Datt=getchandata('batchnotes2','a',[0:1:4],64000,1000);
    %
    pretimems=1000;
    Datt426=getchandata('batchnotes','a',[0:1:4],64000,pretimems);
    preHITms=0;
    fv=findwnoteSPK('batchnotes','a','','',0,[2000 2700],1e4,1,'obs0',0,preHITms);
            clear isTRIG isCATCH
            for i=1:length(fv);isTRIG(i)=fv(i).TRIG;end
            for i=1:length(fv);isCATCH(i)=fv(i).CATCH(1);end
            isHIT1=~isCATCH(find(isTRIG==1)); % for dividing targtimesWN
            isHIT2=find(isTRIG & ~isCATCH); % for identifying fvals
            isESC2=find(~isTRIG | isCATCH);
    
    [avA,t,f]=get_avn('batch422files','a',1,1,'','','obs0'); 
    fs=32000;
    avgwin=250; % points
    clear  rav1 rav2 ravesc ravhit
    currdata=Datt422;
    for i=1:length(currdata);
        i
        rav(i,:)=runningaverage(mean(abs(currdata(i).data(1:53,:))),avgwin);
        %rav2(i,:)=runningaverage(mean(abs(currdata(i).data(54:end,:))),avgwin);        
%         ravesc(i,:)=runningaverage(mean(abs(currdata(i).data(isESC2,:))),avgwin);        
%         ravhit(i,:)=runningaverage(mean(abs(currdata(i).data(isHIT2,:))),avgwin);
    end
    %
    figure;hold on;
%     subplot(length(currdata)+1,1,1)
%     %imagesc(t,f,log(avA));syn;ylim([0 1e4])
%     xlim([-0.2 0.4])
    for i=1:length(currdata)
        subplot(length(currdata)+1,1,i+1);hold on;
%         plot([avgwin/fs+1/fs:1/fs:size(currdata(i).data,2)/fs],ravesc(i,:))
%         plot([avgwin/fs+1/fs:1/fs:size(currdata(i).data,2)/fs],ravhit(i,:),'r')  
        plot([avgwin/fs+1/fs:1/fs:size(currdata(i).data,2)/fs],rav1(i,:),'b')          
%         plot([36000/fs 36000/fs],[min([ravhit(i,:) ravesc(i,:)])*0.9 max([ravhit(i,:) ravesc(i,:)])*1.1],'k')
       plot([36000/fs 36000/fs],[min([rav1(i,:)])*0.9 max([rav1(i,:)])*1.1],'k')        
%         xlim([0.8 1.4]); ylim([min([ravhit(i,:) ravesc(i,:)])*0.9 max([ravhit(i,:) ravesc(i,:)])*1.1])
        xlim([0.8 1.4]); ylim([min([rav1(i,:)])*0.9 max([rav1(i,:)])*1.1])        
    end
%
pitchAll2=pitchcontour(fv,2000,3000);
pts=[-0.8:0.01:0.8];
clear ga
for i=1:length(pts)
    cc=corrcoef(pitchvals,mean(abs(currdata(1).data(:,1.12*fs+pts(i)*fs-0.1*fs:1.12*fs+pts(i)*fs))'));
    ga(i)=cc(2);
end
ind3=find(pitchAll(900,:)>2220 & pitchAll(900,:)<2480);
ind3=ind3(1:120);
pitchvals=mean(pitchAll2(900:1100,:));
figure;plot(pitchAll(900,ind3),mean(abs(Datt422(2).data(ind3,1.12*fs-0.1*fs:1.12*fs))'),'.')
corrcoef(pitchAll(900,ind3),mean(abs(Datt422(2).data(ind3,1.12*fs-0.1*fs:1.12*fs))'))
figure;plot(pitchAll(900,ind3),mean(abs(Datt422(3).data(ind3,1.12*fs-0.1*fs:1.12*fs))'),'.')
corrcoef(pitchAll(900,ind3),mean(abs(Datt422(3).data(ind3,1.12*fs-0.1*fs:1.12*fs))'))
