cd /canary/SY/w62o41/102311_ClickWN
% go into a data directory
% all commands are from matlab

%make a batch file of all .cbin files
!ls -1 *cbin > batch

%pull out  10% of the files wiuievtafsim('batch.catch',templaB);ll make a file called batch.rand
% and batch.notrand (other 90%)
%randsamp('batch',0.1);

% go in and label your target note in batch.rand
% use evsonganaly

%when done labeling -> get average spectrgram around target note

% 'a' is the target note not context notes ('' '')
[avnB,t,f]=get_avn('batch','b',0.2,0.2,'','','obs0'); 

% if 'b' is the target note and you only want 'abc' not 'abd' copntext:
%[avn,t,f]=get_avn('batch.train','b',0.2,0.2,'a','c','obs0'); 
figure
imagesc(t,f,log(avnB));syn;ylim([0,1e4]);

%find time where you want to pull the template out 
[y,i]=min(abs(t-.1));
%the variable 'i' has the columno

%or just plot image without t and freq axis
figure
imagesc(log(avnB));syn;

% pull out template slice  if i chose time slice 94
templaB=mean(avnB(1:2:256,68),2);
size(templaB) %should be 128 point vector  (128 by 1)

%normalize template
templaB(1:6)=0;

templaB=templaB./max(templaB);
%check out template
uievtafsim('batch',templaB);

%get an idea of the counter values and thresholds you want to use


%generate X.tmp files
%   templ is the template vector, 2 is the pre time of the files
%   get the 2 from the rec file use the same # as the T Before value
mk_tempf('batch',templaB,2,'obs0');

%build cntrng struct array
%cntrng(index) -> index is the template # if you have one template it's one
%		 template it only equal to one
%This is set up for case where there are three templates.

% Try one set of cntrng values.


cntrng(1).MIN=2;
cntrng(1).MAX=3;
cntrng(1).NOT=0;
cntrng(1).MODE=1;
cntrng(1).TH=2;
cntrng(1).AND=0;
cntrng(1).BTMIN=0;



cntrng(2).MIN=1;
cntrng(2).MAX=2;
cntrng(2).NOT=0;
cntrng(2).MODE=1;
cntrng(2).TH=2.5;
cntrng(2).AND=1;
cntrng(2).BTMIN=0;



cntrng(3).MIN=1;
cntrng(3).MAX=2;
cntrng(3).NOT=0;
cntrng(3).MODE=1;
cntrng(3).TH=2.5;
cntrng(3).AND=0;
cntrng(3).BTMIN=0


cntrng(4).MIN=2;
cntrng(4).MAX=3;
cntrng(4).NOT=0;
cntrng(4).MODE=1;
cntrng(4).TH=2.2;
cntrng(4).AND=0;
cntrng(4).BTMIN=0

cntrng(5).MIN=2;
cntrng(5).MAX=3;
cntrng(5).NOT=0;
cntrng(5).MODE=1;
cntrng(5).TH=2.2;
cntrng(5).AND=1;
cntrng(5).BTMIN=0




cntrng(1).MIN=3;
cntrng(1).MAX=4;
cntrng(1).NOT=0;
cntrng(1).MODE=1;
cntrng(1).TH=2.2;
cntrng(1).AND=0;
cntrng(1).BTMIN=0



cntrng(2).MIN=3;
cntrng(2).MAX=4;
cntrng(2).NOT=0;
cntrng(2).MODE=1;
cntrng(2).TH=2.4;
cntrng(2).AND=0;
cntrng(2).BTMIN=0



%example if you has a second template
%do a simulation of the counter ranges to see where it would have triggered
get_trigt2('batch',cntrng,0.4,128,1,1); %0.4 = refractory period (sec)


%do a simulation of the coevsunter ranges to see where it would have triggered
get_trigt('batch.catch',cntrng,0.5,128,1,0);

%how well did this template match
[vals,trigs]=triglabel('batch','b',1,1,0,1);

sum(vals)

   232        234        232
% N match    N note    n trig

%More systematic approach.
% another way to determine the threshold
% loop over a bunch of threshold vales pull out hit rates
tvals=[];
for THA=3.5:-0.2:2.5
for TH=3.1:-.1:2.5
	cntrng2(1).TH=THA;
    cntrng2(2).TH=TH;
	get_trigt2('batchT',cntrng2,0.3,128,1,1);
	[vals]=triglabel('batchT','a',1,1,0,1);
    [v,trigs]=triglabel('batchT','a',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
        toff=[toff;trigs(ii).toffset];
    end
	tvals=[tvals;TH,sum(vals),std(toff)];
end

plot(tvals(:,1),tvals(:,2)./tvals(:,3),'bs-')

%pick a TH out of tvals and run triglabel with it
cntrng(1).TH=2.5;
get_trigt('batch2.catch.keep',cntrng,0.3,128,1,1);
[vals,trigs]=triglabel('batchTEMPtest','a',1,1,0,1);

%pull the offsets out of trigs
toff=[];
for ii=1:length(trigs)
	toff=[toff;trigs(ii).toffset];
end
std(toff)
%toff is the trigger offset

figure
imagesc(t,f,log(avna));syn;ylim([0,1e4]);
hold on
plot(hist(toff*1e-3),2000,'k^');
plot((mean(toff)+[-1,1]*std(toff))*1e-3,2000*[1,1],'k-');

figure
imagesc(t,f,log(avna));syn;ylim([0,1e4]);
hold on
plot((toff*1e-3),2000,'*');

%when you decide on a template

wrt_templ('w62o41_temp1025.dat',templaB);

