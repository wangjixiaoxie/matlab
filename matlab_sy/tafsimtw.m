%tafsimtw
%for making templates.

% 'a' is the target note not context notes ('' '')
[avnA,t,f]=get_avn('batch_02.keep.catch','a',0.2,0.2,'','','obs0'); 

%plot with time, frequency
figure
imagesc(t,f,log(avnA));syn;ylim([0,1e4]);

%for template making I plot image without t and freq axis
figure
imagesc(log(avnA));syn;

% I choose slice to pull template from.
templa1=mean(avnA(1:2:256,64:65),2);
size(templa1) %should be 128 point vector  (128 by 1)

%normalize first six values of template to 0.
templa1(1:6)=0;
templa1=templa1./max(templa1);
%check out template with gui
uievtafsim('batch',templa1);

%get an idea of the counter values and thresholds you want to use


%generate X.tmp files
%   templ is the template vector, 2 is the pre time of the files
%   get the 2 from the rec file use the same # as the T Before value
mk_tempf('batch_02.keep.catch',bk90bk46,2,'obs0');

%build cntrng struct array
%cntrng(index) -> index is the template # if you have one template it's one
%		 template it only equal to one
%This is set up for case where there are three templates.

% Try one set of cntrng values.

%min threshold
cntrng(1).MIN=1;
cntrng(1).MAX=2;
%true/false logic, true->note=0
cntrng(1).NOT=0;
%evtafmode=1; birdtafmode=0;
cntrng(1).MODE=1;
%threshold
cntrng(1).TH=2.2;
%and/or logic with other templates.
cntrng(1).AND=0;
cntrng(1).BTMIN=0



cntrng(2).MIN=1;
cntrng(2).MAX=2;
cntrng(2).NOT=0;
cntrng(2).MODE=1;
cntrng(2).TH=2.2;
cntrng(2).AND=0;
cntrng(2).BTMIN=0

% 
% 
cntrng(3).MIN=1;
cntrng(3).MAX=2;
cntrng(3).NOT=0;
cntrng(3).MODE=1;
cntrng(3).TH=2.2;
cntrng(3).AND=1;
cntrng(3).BTMIN=0
% 
% 
cntrng(4).MIN=2;
cntrng(4).MAX=15;
cntrng(4).NOT=0;
cntrng(4).MODE=0;
cntrng(4).TH=1.9;
cntrng(4).AND=0;
cntrng(4).BTMIN=0
% 
cntrng(5).MIN=18;
cntrng(5).MAX=32;
cntrng(5).NOT=0;
cntrng(5).MODE=0;
cntrng(5).TH=3;
cntrng(5).AND=1;
cntrng(5).BTMIN=0
% 
% cntrng(6).MIN=3;
% cntrng(6).MAX=4;
% cntrng(6).NOT=0;
% cntrng(6).MODE=1;
% cntrng(6).TH=2.6;
% cntrng(6).AND=0;
% cntrng(6).BTMIN=0
% 
% cntrng(7).MIN=3;
% cntrng(7).MAX=4;
% cntrng(7).NOT=0;
% cntrng(7).MODE=1;
% cntrng(7).TH=2.6;
% cntrng(7).AND=1;
% cntrng(7).BTMIN=0
% 
% 

%SIMULATION OF TRIGGERING, with those cntrng values.
%example if you has a second template
%do a simulation of the counter ranges to see where it would have triggered
get_trigt2('batch_02.keep.catch',cntrng,0.085,128,1,1);


%how well did this template match (addx 1 for simulation, 0 for real)
[vals,trigs]=triglabel('batch_02.keep.catch','b',1,1,0,1);

sum(vals)

   232        234        232
% N match    N note    n trig

%More systematic approach.
% another way to determine the threshold
% loop over a bunch of threshold vales pull out hit rates
tvals=[];
for TH=2.8:-.1:2.5[avnA,t,f]=get_avn('batch_02.keep.catch','a',0.2,0.2,'','','obs0'); 
	cntrng(1).TH=TH;
    cntrng(2).TH=TH;
    cntrng(3).TH=TH
	get_trigt2('batch04.keep',cntrng,0.18,128,1,1);
	[vals]=triglabel('batch04.keep','b',1,1,0,1);
	tvals=[tvals;TH,sum(vals)];
 end

figure, plot(tvals(:,1),tvals(:,2)./tvals(:,3),'bs-')

%pick a TH out of tvals and run triglabel with it
cntrng(1).TH=2;
get_trigt('batch_02.keep.catch',cntrng,0.085,128,1,1);  
[vals,trigs]=triglabel('batch_02.keep.catch','b',1,1,0,1);

%pull the offsets out of trigs
toff=[];
for ii=1:length(trigs)
	toff=[toff;trigs(ii).toffset];
end

%toff is the trigger offset

figure
imagesc(t,f,log(avna));syn;ylim([0,1e4]);
hold on
plot(mean(toff*1e-3),2000,'k^');
plot((mean(toff)+[-1,1]*std(toff))*1e-3,2000*[1,1],'k-');


%THIS CODE WRITES THE TEMPLATE TO A DATA FILE

wrt_templ('bk90bk46.dat',tmpcomb);
