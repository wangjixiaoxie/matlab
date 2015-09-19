% Get average spectrgram around target note
randsamp('batch',0.1)
        % 'a' is the target note not context notes ('' '')
        [avA,t,f]=get_avn('batch','a',1,1,'','','obs0'); 
                % if 'b' is the target note and you only want 'abc' not 'abd' copntext:
                %[avn,t,f]=get_avn('batch.train','b',0.2,0.2,'a','c','obs0'); 
        figure
        imagesc(t,f,log(avA));syn;ylim([0,1e4]);xlim([-1 1])

%% TEMPLATE %%%
% Look at it without time axis to help make template
    figure
    imagesc(log(avA));syn;
% pull out template slice  if i chose time slice 67
templa1=mean(avA(1:2:256,[70]),2);
size(templa1) %should be 128 point vector  (128 by 1)
%normalize template
templa1(1:6)=0;
templa1=templa1./max(templa1);
%check out template
uievtafsim('batch.rand',templaStim);

%get an idea of the counter values and thresholds you want to use
%generate X.tmp files
%   templ is the template vector, 2 is the pre time of the files
%   get the 2 from the rec file use the same # as the T Before value
mk_tempf('batch.rand',templaStim,2,'obs0');

%build cntrng struct array
%cntrng(index) -> index is the template # if you have one template it's one
%		 template it only equal to one
%This is set up for case where there are three templates.

% Try one set of cntrng values.

    cntrng1(1).MIN=5;
    cntrng1(1).MAX=15;
    cntrng1(1).NOT=0;
    cntrng1(1).MODE=0;
    cntrng1(1).TH=2.6;
    cntrng1(1).AND=1;
    cntrng1(1).BTMIN=0;

  % example if you have a second template
    cntrng1(2).MIN=2;
    cntrng1(2).MAX=3;
    cntrng1(2).NOT=0;
    cntrng1(2).MODE=1;
    cntrng1(2).TH=2;
    cntrng1(2).AND=1;
    cntrng1(2).BTMIN=0;


%do a simulation of the counter ranges to see where it would have triggered
get_trigt2('batch.rand',cntrng1,0.2,128,1,1);

%how well did this template match
[vals,trigs]=triglabel('batch.rand','b',1,1,0,1);

sum(vals)

%   232        234        232
% N match    N note    n trig

%pull the offsets out of trigs
toff=[];
for ii=1:length(trigs)
	toff=[toff;trigs(ii).toffset];
end
std(toff)  %toff is the trigger offset


% More methodolical approach

tvals=[];

for TH=3:-.2:2
    cntrng1(1).TH=TH;
    get_trigt2('batch.rand',cntrng1,0.2,128,1,1);
	[vals]=triglabel('batch.rand','b',1,1,0,1);
    [v,trigs]=triglabel('batch.rand','b',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
        toff=[toff;trigs(ii).toffset];
    end
	tvals=[tvals;TH,sum(vals),std(toff)];
end
tvals

%when you decide on a template

  wrt_templ('bk48w74stim_72710.dat',templaStim);

