% go into a data directory
% all commands are from matlab

%make a batch file of all .cbin files
!ls -1 *cbin > batch
dirf('*.cbin','batch720Afiles');
%pull out  10% of the files will make a file called batch.rand
% and batch.notrand (other 90%)
randsamp('batch',0.1);
% Make a batchNOTES file by find and replace all .cbin with .cbin.not.mat

% go in and label your target note in batch.rand
% use evsonganaly

%when done labeling -> get average spectrgram around target note

%batch.rand = batchTEMPtest in this example
% 'a' is the target note not context notes ('' '')
[avZ,t,f]=get_avn('batch.catch.rand','z',0.2,0.2,'','','obs0'); 

% if 'b' is the target note and you only want 'abc' not 'abd' copntext:
%[avn,t,f]=get_avn('batch.train','b',0.2,0.2,'a','c','obs0'); 
figure
imagesc(t,f,log(avC));syn;ylim([0,1e4]);xlim([-0.05 0.45])

%find time where you want to pull the template out 
[y,i]=min(abs(t-.1));
%the variable 'i' has the columno

%or just plot image without t and freq axis
figure
imagesc(log(avZ));syn;
%%TEMPLATE%%
% pull out template slice  if i chose time slice 94
templa5=mean(avZ(1:2:256,[64]),2);
size(templa5) %should be 128 point vector  (128 by 1)
%normalize template
templa5(1:6)=0;
templa5=templa5./max(templa5);
%check out template
uievtafsim('batch.catch.rand',templa5);

%get an idea of the counter values and thresholds you want to use


%generate X.tmp files
%   templ is the template vector, 2 is the pre time of the files
%   get the 2 from the rec file use the same # as the T Before value
mk_tempf('batch.catch.rand',templaD,2,'obs0');
pitch1=jc_PitchData208(shifted,1024,1020,2,3000,4200,1,'obs0');
%build cntrng struct array
%cntrng(index) -> index is the template # if you have one template it's one
%		 template it only equal to one
%This is set up for case where there are three templates.

% Try one set of cntrng values.

cntrngA(1).MIN=2;
cntrngA(1).MAX=3;
cntrngA(1).NOT=0;
cntrngA(1).MODE=1;
cntrngA(1).TH=2.5;
cntrngA(1).AND=1;
cntrngA(1).BTMIN=0;


cntrngA(2).MIN=1;
cntrngA(2).MAX=10;
cntrngA(2).NOT=1;
cntrngA(2).MODE=1;
cntrngA(2).TH=2;
cntrngA(2).AND=1;
cntrngA(2).BTMIN=0;


cntrng4(3).MIN=4;
cntrng4(3).MAX=1down;
cntrng4(3).NOT=0;
cntrng4(3).MODE=1;
cntrng4(3).TH=3.6;
cntrng4(3).AND=0;
cntrng4(3).BTMIN=0;
%example if you has a second template
%do a simulation of the counter ranges to see where it would have triggered
get_trigt2('batch.catch.rand',cntrngD,0.4,128,1,1);

%do a simulation of the coevsunter ranges to see where it would have triggered
%get_trigt('batchJC',cntrng,0.1down,128,1,0);

%how well did this template match
[vals,trigs]=triglabel('batch.catch.rand','d',1,1,0,1);

sum(vals)

%   232        234        232
% N match    N note    n trig

%More systematic approaches
% another way to determine the threshold
% loop over a bunch of threshold vales pull out hit rates
% and standevs - no need to plot
%%% vary the threshold
suplot(tvals(:,1),tvals(:,2)./tvals(:,3),'bs-')
%jc_IOtrigger shows how to vary another parameter - such as refract period
tvals=jc_IOtrigger(cntrng);

% Try one set of cntrng values.

cntrngA(1).MIN=1;
cntrngA(1).MAX=2;
cntrngA(1).NOT=0;
cntrngA(1).MODE=0;
cntrngA(1).TH=2;
cntrngA(1).AND=1;
cntrngA(1).BTMIN=0;
%pick a TH out of tvals and run triglabel with it
cntrng3(1).TH=2.3;
get_trigt('batchT',cntrng3,0.3,128,1,1);
[vals,trigs]=triglabel('batchTMP2files','a',1,1,0,1);

%pull the offsets out of trigs
toff=[];
for ii=1:length(trigs)
	toff=[toff;trigs(ii).toffset];
end
std(toff)  %toff is the trigger offset

%PLOT offsets against the spectrogram of the note
figure
imagesc(t,f,log(avna));syn;ylim([0,1e4]);
hold on
plot((toff*1e-3),2000,'*');
[vals,trigs]=triglabel('batchJCfiles','a',1,1,0,1);

sum(vals)

% Shows you the real std when you factor in alignment with jc_AlignCT
% Also gives you the fvalsfor the next steps
Contingency(toff,batchNOTE) 

tvals=[];

for TH=3:-.1:2

    cntrngA(1).TH=TH;
    get_trigt2('batchJCfiles',cntrngA,0.2,128,1,1);
	[vals]=triglabel('batchJCfiles','a',1,1,0,1);
    [v,trigs]=triglabel('batchJCfiles','a',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
        toff=[toff;trigs(ii).toffset];
    end
	tvals=[tvals;TH,sum(vals),std(toff)];
end
tvals

%What is the expected shift based on variation in the data
shifted=jc_AlignCT(fvals,toff); %toff
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,191down0,2600,[1 2 3],'obs0',1);
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

  wrt_templ('templaD_1031.dat',templaD);

