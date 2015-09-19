% Control #1 - How does learning reflect what was actually reinforced as
% opposed to what you might expect to be reinforced?

% Reviewer #2 suggested an experiment in which aversive reinforcement would be
% delivered not only to all performances below a threshold criterion 
% (as in the original experiments) but also to a random subset of
% performances above that criterion.  If learning reflected the history of
% associations between performance and reinforcement, as we proposed, then
% the random subset of performances above the criterion that escaped
% aversive reinforcement should predict learning better than the random 
% subset of performances above the criterion that received aversive reinforcement.
% In constrast, if learning only reflects parameters of the motor control apparatus
% and does not reflect the specific history of performance and
% reinforcement, then the performances that randomly received reinforcement
% (random escapes) should predict learning better than the performances that randomly
% received reinforcement (random hits).

% To determine whether learning reflects the history of recent associations
% between performance and reinforcement, we employed a trial-by-trial
% analysis of learning during the reinforcement/training period.  
% Since white noise playback reduces our ability to quantify fundamental 
% frequency performance after the contingency time, we specifically analyzed 
% performance before (and during) the contingency time.
% On each trial (i.e. performance of a syllable during the reinforcement period), 
% learning was defined as difference from mean performance at baseline
% (before any reinforcement).  We performed linear regression to compare
% how well the fundamental frequency performance in a given trial (learning) could be
% predicted by fundamental frequency performance on recent random hits and
% recent random escapes.  If learning reflects the history of associations 
% between performance and reinforcement, as we proposed, then recent random
% escapes should predict learning better than random recent hits.

% Random recent escapes predicted learning better than random recent hits.
% We performed n=4 experiments in 2 birds.  In each experiment, we analyzed 
% at least one day of song, containing at least 900 "trials."  Results were
% consistent across experiments: in each recent time window, random escapes
% predicted learning better than random hits.

load /cardinal3/TrialbyTrial_Sept2010.mat


window2=[901:1101];
window5=[301:501];
whichnotes2=[1:size(Experiment(2).pitchWNon,2)-500];   % 1:size(Experiment(2).pitchWNon,2);
whichnotes5=[1:size(Experiment(5).pitchWNon,2)-500];   % 1:size(Experiment(5).pitchWNon,2);
Pred2=(mean(Experiment(2).pitchWNon(window2,min(whichnotes2)-1+find(Experiment(2).indEscapeAbove(whichnotes2)))'))-mean(Experiment(2).pitchBaseline(window2,:)');
PredH2=(mean(Experiment(2).pitchWNon(window2,min(whichnotes2)-1+find(Experiment(2).indHitAbove(whichnotes2)))'))-mean(Experiment(2).pitchBaseline(window2,:)');
Act2=(mean(pitchPost2(window2-600,:)')-mean(Experiment(2).pitchBaseline(window2,:)'));
Pred5=(mean(Experiment(5).pitchWNon(window5,min(whichnotes2)-1+find(Experiment(5).indEscapeAbove(whichnotes5)))'))-mean(Experiment(5).pitchBaseline(window5,:)');
PredH5=(mean(Experiment(5).pitchWNon(window5,min(whichnotes2)-1+find(Experiment(5).indHitAbove(whichnotes5)))'))-mean(Experiment(5).pitchBaseline(window5,:)');
Act5=mean(pitchPost5(window5,:)')-mean(Experiment(5).pitchBaseline(window5,:)');


% Column 1 - all --- Column 2 - ignore most recent 500
corrcoef(Pred5,Act5) % 0.73  % 0.71
corrcoef(Pred2,Act2) % 0.89  % 0.90
corrcoef(Pred2,Act5) % 0.55  % 0.51
corrcoef(Pred5,Act2) % 0.70  % 0.56

%%% Another bird
window6=[350:550];
whichnotes6= 1:size(Experiment(6).pitchWNon,2);
Pred6=(mean(Experiment(6).pitchWNon(window6,min(whichnotes6)-1+find(Experiment(6).indEscapeAbove(whichnotes6)))'))-mean(Experiment(6).pitchBaseline(window6,:)');
PredH6=(mean(Experiment(6).pitchWNon(window6,min(whichnotes6)-1+find(Experiment(6).indHitAbove(whichnotes6)))'))-mean(Experiment(6).pitchBaseline(window6,:)');
Act6=(mean(pitchPost6(window6,:)')-mean(Experiment(6).pitchBaseline(window6,:)'));
corrcoef(Pred6,Act6) % 0.8469

figure;hold on;
subplot(221);hold on;
plot(Pred2,'b');plot(Act2,'r')
subplot(222);hold on;
plot(Pred5,'b');plot(Act5,'r')
subplot(223);hold on;
plot(Pred6,'b');plot(Act6,'r')

plot(mean(pitchPost2'))
plot(mean(pitchPost5'),'r')
plot(mean(Experiment(5).pitchWNon(:,find(Experiment(5).indEscapeAbove))'),'r')
plot(mean(Experiment(2).pitchWNon(:,find(Experiment(2).indEscapeAbove))'),'b')
%%%%%%%%
%%%%%%%%
%%% Code 
%%%%%%%%
%%%%%%%%
% Make sure auto-segmentation is okay
    evsonganaly
% autolabel
    mk_tempf('batch922.rand',templaB,2,'obs0');
    get_trigt2('batch922.rand',cntrngB,0.2,128,1,1);
    note='b';
    label_trigs('batch',note,'obs0',10000,1,1,5,30);
% Check autolabels
    evsonganaly
% Make a 'batchnotes' file
% Record information about experiment
    expnumber=7; % Which experiment overall?
    Experiment(expnumber).bird='r87g80';
    Experiment(expnumber).expnum=2; % Which experiment for this bird?
    Experiment(expnumber).wnday='Sept_24_2010'; % First day of WN (8:30pm doesn't count)
    Experiment(expnumber).note=note;
% Process data
    Experiment(expnumber).targeting=gettarg('batch924files',Experiment(expnumber).note);
    fvBaseline=findwnoteJC('batchnotes','b','a','',0,[6000 8100],8500,1,'obs0',0);
    fvWNon=findwnoteNN('batch924notes','b','','',0,[6000 8100],8500,1,'obs0',0);
    [Experiment]=trialbytrialpre(Experiment,expnumber,fvBaseline,fvWNon);

    
    
%%%% Look at next    
  for i=1:6
    avgposthit(i).data=median(Experiment(i).pitchWNon(:,find(Experiment(i).indHitAbove(1:end-1)==1)+1)');
    avgpostescape(i).data=median(Experiment(i).pitchWNon(:,find(Experiment(i).indEscapeAbove(1:end-1)==1)+1)');
end
figure;hold on;
for i=1:6
    plot(avgpostescape(i).data(Experiment(i).ptwindow)-avgposthit(i).data(Experiment(i).ptwindow))
    m(i)=median(avgpostescape(i).data(Experiment(i).ptwindow)-avgposthit(i).data(Experiment(i).ptwindow));
end
[h,p]=ttest(m)  

%%%%% Doesn't seem promising
for i=1:6
    avghit(i).data=mean(Experiment(i).pitchWNon(:,Experiment(i).indHitAbove(1:200))');
    avgescape(i).data=mean(Experiment(i).pitchWNon(:,Experiment(i).indEscapeAbove(1:200))');
end
figure;hold on;
for i=1:6
    plot(mean(Experiment(i).pitchWNon(Experiment(i).ptwindow,end-100:end)')-avgescape(i).data(Experiment(i).ptwindow))
    plot(mean(Experiment(i).pitchWNon(Experiment(i).ptwindow,end-100:end)')-avghit(i).data(Experiment(i).ptwindow),'r')
end
[h,p]=ttest(m)  

    
%%%%%    
    
    
    
    
% Set parameters for analysis - load processed data here
figure;plot(std(Experiment(expnumber).pitchWNon(:,Experiment(expnumber).indEscapeAbove)'),'b')
hold on;plot(std(Experiment(expnumber).pitchWNon(:,Experiment(expnumber).indHitAbove)'),'r')
    Experiment(expnumber).ptwindow=[260:460]; % should be largest window with similar stan dev for both
    startingptwindow=round(median(Experiment(expnumber).targeting));
    if startingptwindow-64>min(Experiment(expnumber).ptwindow) & startingptwindow<max(Experiment(expnumber).ptwindow)
        Experiment(expnumber).targwindow=[startingptwindow-64:startingptwindow]; % neither is problematic
    else if startingptwindow-64<min(Experiment(expnumber).ptwindow) & startingptwindow>max(Experiment(expnumber).ptwindow)
            Experiment(expnumber).targwindow=Experiment(expnumber).ptwindow; % both are problematic
        else if startingptwindow-64<min(Experiment(expnumber).ptwindow) % small one is problematic
                Experiment(expnumber).targwindow=[min(Experiment(expnumber).ptwindow):startingptwindow];
        else Experiment(expnumber).targwindow=[startingptwindow-64:max(Experiment(expnumber).ptwindow)];
        end
        end
    end
    histwindow=100;
% Analyze data
    % targwindow vs pt window
    % subtraction vs division vs dotproduct
        % subtraction tells us
        % division tells us 
        % dotproduct tells us linear relationship (median of corr coef)
        figure;hold on;
        runwin=20;
%   for expnumber=1:3
for expnumber=1:6
    [fmrecentescapes,fmrecenthits,fmrecent]=trialbytrialcore(Experiment,expnumber,histwindow,'ptwindow','subtraction');
    Experiment(expnumber).escapeSubtraction=fmrecentescapes;
    Experiment(expnumber).hitSubtraction=fmrecenthits;  
%     subplot(1,3,expnumber);hold on;
end
   figure;hold on; plot(runningaverage(Experiment(expnumber).escapeSubtraction,runwin),'b')
    plot(runningaverage(Experiment(expnumber).hitSubtraction,runwin),'r')
% 
%   end

% FOR EXPERIMENT 5, only look at 50-350 (learning w/o noisiness).
clear a
% Is the hitdistance greater than the escape distance?
for i=1:length(Experiment)
  a(i,:)=Experiment(i).hitSubtraction-Experiment(i).escapeSubtraction;
end
% Is this effect less present at the targ window (i.e. does it indicate'shape')?
a=a([1 2 4 5 6],:);
figure;plot(runningaverage(mean(a),20))
hold on;plot(runningaverage(a,20)-runningstd(a,20)/sqrt(length(Experiment)))
  
%     [fmrecentescapes,fmrecenthits,fmrecent]=trialbytrialcore(Experiment,expnumber,histwindow,'ptwindow','division');
%     Experiment(expnumber).escapeDivision=fmrecentescapes;
%     Experiment(expnumber).hitDivision=fmrecenthits;
%     [fmrecentescapes,fmrecenthits,fmrecent]=trialbytrialcore(Experiment,expnumber,histwindow,'ptwindow','corrcoef');    
%     Experiment(expnumber).escapeCorrcoef=fmrecentescapes;
%     Experiment(expnumber).hitCorrcoef=fmrecenthits;   
  
    
% PLOT
expnumber=2
expnumber=3
expnumber=3
figure;hold on;
runwin=20;
    subplot(133);hold on;
    plot(runningaverage(Experiment(expnumber).escapeSubtraction,runwin),'b')
    plot(runningaverage(Experiment(expnumber).hitSubtraction,runwin),'r')
%     subplot(132);hold on;
%     plot(runningaverage(Experiment(expnumber).escapeDivision,runwin),'b')
%     plot(runningaverage(Experiment(expnumber).hitDivision,runwin),'r')
%     subplot(133);hold on;
%     plot(runningaverage(Experiment(expnumber).escapeCorrcoef,runwin),'b')
%     plot(runningaverage(Experiment(expnumber).hitCorrcoef,runwin),'r')
% 

% bk80w28 - data recorded on launchpad
% Experiment 1
    % baseline: /cardinal3/bk80w28/0903_baseline/batchnotes
        % a random sampling of songs from 9.07.10
    % wn: /cardinal3/bk80w28/0907_wnon/batchnotes
        % every song in the folder - starts at dark on 9.07.10 and ends at 10am on 9.09.10
% Experiment 2
    % baseline: /cardinal3/bk80w28/0910_wnoff/batchnotes
        % a random sampling of songs from 9.12.10
    % wn: /cardinal3/bk80w28/0913_wnon/batchnotes
        % every song in the folder - starts at dark on 9.12.10 and ends at 4:12pm on 9.13.10

% bk91w60 - data recorded on bigbird? - computer #4 on Evren's desk
    
    % Experiment #1
    load /cardinal3/bk91w60/NatNeuroCTL091710.mat
    edit bk91w60ctls.m
    % window is [260:310] - before WN
    % First point is most recent.

    
    
% Control #2 - How profound are the pitch trajectory changes that can be
% elicited by targeted WN?  Can we turn a stack note into a sweep? 

% Look for great examples...