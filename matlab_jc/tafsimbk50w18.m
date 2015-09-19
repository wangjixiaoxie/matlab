% go into a data directory
% all commands are from matlab

%make a batch file of all .cbin files
!ls -1 *cbin > batch
dirf('*.cbin','batch');
%pull out  10% of the files will make a file called batch.rand
% and batch.notrand (other 90%)
randsamp('batch',0.1);
% Make a batchNOTES file by find and replace all .cbin with .cbin.not.mat

% go in and label your target note in batch.rand
% use evsonganaly

%when done labeling -> get average spectrgram around target note

%batch.rand = batchTEMPtest in this example
% 'a' is the target note not context notes ('' '')
[avB,tb,fb]=get_avn('batchtest','a',0.2,0.2,'','','obs0'); 

% if 'b' is the target note and you only want 'abc' not 'abd' copntext:
%[avn,t,f]=get_avn('batch.train','b',0.2,0.2,'a','c','obs0'); 
figure
imagesc(tb,fb,log(avB));syn;ylim([0,1e4]);

%find time where you want to pull the template out 
[y,i]=min(abs(t-.1));
%the variable 'i' has the columno

%or just plot image without t and freq axis
figure
imagesc(log(avna));syn;
%%TEMPLATE%%
% pull out template slice  if i chose time slice 94
templa=mean(avB(1:2:256,[90]),2);
size(templa) %should be 128 point vector  (128 by 1)
%normalize template
templa(1:6)=0;
templa=templa./max(templa);
%check out template
uievtafsim('batchtest',templa);

%get an idea of the counter values and thresholds you want to use


%generate X.tmp files
%   templ is the template vector, 2 is the pre time of the files
%   get the 2 from the rec file use the same # as the T Before value
mk_tempf('batchtest',templa,2,'obs0');

%build cntrng struct array
%cntrng(index) -> index is the template # if you have one template it's one
%		 template it only equal to one
%This is set up for case where there are three templates.

% Try one set of cntrng values.

cntrng1(1).MIN=5;
cntrng1(1).MAX=6;
cntrng1(1).NOT=0;
cntrng1(1).MODE=0;
cntrng1(1).TH=4;
cntrng1(1).AND=0;
cntrng1(1).BTMIN=0;


cntrng(2).MIN=6;
cntrng(2).MAX=7;
cntrng(2).NOT=0;
cntrng(2).MODE=0;
cntrng(2).TH=4.4;
cntrng(2).AND=0;
cntrng(2).BTMIN=0;

%example if you has a second template
%do a simulation of the counter ranges to see where it would have triggered
get_trigt2('batchtest',cntrng1,0.05,128,1,1);

%do a simulation of the coevsunter ranges to see where it would have triggered
%get_trigt('batchJC',cntrng,0.5,128,1,0);

%how well did this template match
[vals,trigs]=triglabel('batchtest','a',1,1,1,1);

sum(vals)

   232        234        232
% N match    N note    n trig

%More systematic approaches
% another way to determine the threshold
% loop over a bunch of threshold vales pull out hit rates
% and standevs - no need to plot
%%% vary the threshold
tvals=[];
for TH=5:-.2:3.4
	cntrng(1).TH=TH;
	get_trigt2('batchtest',cntrng,0.3,128,1,1);
	[vals,trigs]=triglabel('batchtest','a',1,1,1,1);
    toff=[];
    for ii=1:length(trigs)
        toff=[toff;trigs(ii).toffset];
    end
    tvals=[tvals;TH,sum(vals),std(toff)];
end
plot(tvals(:,1),tvals(:,2)./tvals(:,3),'bs-')
%jc_IOtrigger shows how to vary another parameter - such as refract period
tvals=jc_IOtrigger(cntrng);


%pick a TH out of tvals and run triglabel with it
cntrng(1).TH=4.4;
get_trigt('batchtest',cntrng,0.05,128,1,1);
[vals,trigs]=triglabel('batchtest','a',1,1,1,1);

%pull the offsets out of trigs
toff=[];
for ii=1:length(trigs)
	toff=[toff;trigs(ii).toffset];
end
std(toff)  %toff is the trigger offset

%PLOT offsets against the spectrogram of the note
figure
imagesc(tb,fb,log(avB));syn;ylim([0,1e4]);
hold on
plot((toff*1e-3),2000,'*');

% Shows you the real std when you factor in alignment with jc_AlignCT
% Also gives you the fvalsfor the next steps
Contingency(toff,batchNOTE) 


%What is the expected shift based on variation in the data
shifted=jc_AlignCT(fvals,toff); %toff
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,1950,2600,[1 2 3],'obs0',1);
target=jc_contingencyAV1(pitch,'obs0',1,'above',[430 497]);
figure;plot(avg);hold on;plot(target,'r')

%Do the results from my pitch estimate generally agree with those from
%the online fft estimator?
pitchPOS=pitch(:,indX); %ignore false negatives
Bar=2380; %mean(valsX(:,2))+1*std(valsX(:,2));
best_guess=mean(pitchPOS(430:497,:));
figure;hold on; 
plot(best_guess(:),1,'*','Color','k')
plot(best_guess(ind),1,'*','Color','r')


%when you decide on a template

wrt_templ('bk50w18_092808.dat',templa);

