% go into a data directory
% all commands are from matlab


%STEP 1: LABEL
%STEP 2: summary_stats2.m
[fvalsT,FFstatsT,avgT,toffsT]=summary_stats2(templa,'batchtest','batchtestnotes',cntrng);
%STEP 3: summary_plot.m
[xax]=summary_plot(FFs,batch,avgpitch);
%STEP 4: plot a baseline curve
hold on;plot(xax,a913,'k')





%make a batch file of all .cbin file
dirf('*.cbin','batch');
%pull out  10% of the files will make a file called batch.rand
% and batch.notrand (other 90%)
randsamp('batch',0.1);
% Make a batchNOTES file by find and replace all .cbin with .cbin.not.mat

% go in and label your target note in batch.rand
evsonganaly

%%TEMPLATE%%
%load file containing the template
% Load file containing cntrng

%generate X.tmp files:
%   templa is the template vector, 2 is the pre time of the files
%   get the 2 from the rec file use the same # as the T Before value
mk_tempf('batchJC',templa,2,'obs0');

%generate X.rec files:
%do a simulation of the counter ranges to see where it would have triggered
get_trigt2('batchJC',cntrng,0.3,128,1,1);

%%%fvals with real triggers 
fvals=findwnote4('batchJCnotes','a','','',0,[2000 2700],8500,1,'obs0',0);
%%%fvals with triggers from X.rec files
fvalsX=findwnote4('batchJCnotes','a','','',0,[2000 2700],8500,1,'obs0',1);

%how well did this template match - tells you false pos/neg due to
%flaw in template as opposed to contingency
[vals,trigs]=triglabel('batchJC','a',1,1,0,1);
sum(vals)
%   68              72                     77
%Matches        %labeled notes       %Total positives
                                     % Column 1+false positives
toff=[];
for ii=1:length(trigs)
	toff=[toff;trigs(ii).toffset];
end

%Select notes that are not false negatives due to template flaw
%Find indices of these notes
vvalsX=getvals(fvalsX,1,'TRIG'); %Looks in X.rec
indX=find(vvalsX(:,3)~=0);
%Select notes that are negatives due to escaping the contingency
%This will be all the notes that were not false negatives
vvals=getvals(fvals,1,'TRIG'); %Looks in .rec
vv=vvals(indX,:);
ind=find(vv(:,3)==0);
Evals=evtaf_freq('batchJC',[1950 2600],'a',128,'obs0',1,1);
figure; hold on;
%All of these notes should match since this is a simulation 
%of the online contingency
plot(Evals(:,2),'*');plot(ind,Evals(ind,2),'*','Color','r')
%How well is the online frequency calculator actually doing?
%Compare to getvals (more rigorous pitch calculator)
figure;hold on;
plot(vv(:,2),'*')
plot(ind,vv(ind,2),'*','Color','r')


%Plot expected shift based on baseline data.
shifted=jc_AlignCT(fvals,toff); %toff
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,1950,2600,[1 2 3],'obs0',1);
target=jc_contingencyAV1(pitch,'obs0',1,'above',[1059 1122]);
figure;plot(avg);hold on;plot(target,'r')

%Do the results from my pitch estimate generally agree with those from
%the online fft estimator?
pitchPOS=pitch(:,indX); %ignore false negatives
Bar=2380; %mean(valsX(:,2))+1*std(valsX(:,2));
best_guess=mean(pitchPOS(1059:1122,:));
figure;hold on; 
plot(best_guess(:),1,'*','Color','k')
plot(best_guess(ind),1,'*','Color','r')


%Generate white noise file
%unfilteredwn.m explains how to make white noise without a notch
%lowpasswn.m is tim's program for making a notch filter
%read back in a data file so you can pad with zeros
    [data,fs]=wavread('wn50.wav'); 

%%% In depth analysis of variability
normalized=residuals(pitch);
%chop down to size
[xax,totals]=WaveletAnalysisJC2(norm_resids_chopped);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%September 16,2008
findcatch('batch')  % finds the catch trials (the details of which
                    % are written into the rec files