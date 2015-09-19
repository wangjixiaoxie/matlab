%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls


avls.baspath='/oriole/bk63w43/'
avls.datdir='contoursum'

wnin(1).tmon{1}={'2009-8-22 7' '2009-8-24 7'}
wnin(1).tmon{2}={'2009-9-02 7' '2009-9-04 7'}
wnin(1).tmon{3}={'2009-9-26 7'}
wnin(1).tmon{4}={'2009-11-01 7'}
wnin(1).tmoff{1}={'2009-8-24 7' '2009-8-29 7'}
wnin(1).tmoff{2}={'2009-9-04 7' '2009-9-09 7'} 
wnin(1).tmoff{3}={'2009-10-12 7'}
wnin(1).tmoff{4}={'2009-11-10 7'}


wnrevin(1).tmon{1}={'2009-8-29 7'}
wnrevin(1).tmoff{1}={'2009-8-31 7'}
wnrevin(1).tmon{2}={'2009-9-09 7'}
wnrevin(1).tmoff{2}={'2009-9-12 7'}
wnrevin(1).tmon{2}={'2009-9-09 7'}
wnrevin(1).tmoff{2}={'2009-9-12 7'}
wnrevin(1).tmon{3}={'2009-10-12 7'}
wnrevin(1).tmoff{3}={'2009-10-15 7'}
wnrevin(1).tmon{4}={'2009-11-11 7'}
wnrevin(1).tmoff{4}={'2009-11-15 7'}


avls.wnin=wnin
avls.wnrevin=wnrevin

avls.NT{1}='a';
avls.NT{2}='a';
avls.NT{3}='a';
pathvl{1}='ac1212/'
catchvl{1}='batch13.keepb';tmn{25}=['7'];tmf{25}=['10:30'];

pathvl{2}='2001213/'
catchvl{2}='batch.keepb';tmn{26}=['10:30'];tmf{26}=['13:08'];

pathvl{3}='ac1213/'
catchvl{3}='batchb';tmn{27}=['13:08'];tmf{27}=['21'];
% avls.pvls{2}='400mu1204/'
% avls.cvl{2}='batch.keep';
% avls.NT{2}='a'
% 
% avls.pvls{3}='ac1204/'
% avls.cvl{3}='batch04comb';
% avls.NT{3}='a'
avls.pvls=pathvl;
avls.cvl=catchvl;
for ii=1:8
    
    avls.del{ii}=.065
end
for ii=9
    avls.del{ii}=.095;
end
%do pitch analysis on these batch inds.

avls.mkfv=[3]

    avls.con_tbinshft=-0.2;
    %shift for pitch analysis
    avls.pt_tbinshft=.037;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  %shift for contours
   
    avls.fbins=[2000 3000]
%     avls.conbins=[2400 2800];
    
    
    avls.basruns=[]
    avls.contanal=1
avls.contfv=[1:3]
avls.ptfv=[1:3]
avls.analfv=[3]
avls.catchstimfv=[]
avls.contnum=1:3;
    avls.con_tbinshft=-.01;
    %shift for pitch analysis
    avls.pt_tbinshft=.079;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  %shift for contours
    avls.mnbas=[7172]
    avls.stdbas=[118]
    avls.basruns=1;
    
    avls.SOUNDSTR='wav'
    avls.STIMSTR='rig'
    avls.STIMBNDS=[200 0]
     avls.HST_EDGES=[6300:50:8000]
%     avls.STAN_RUNS=[1:7 13:15 17:19 23:25 27:36 39 41:52 54:56]
    avls.STAN_RUNS=[1:3]
    conbinsmid=[30]	
    conbinshi=[73:87]
    conbinsnorm=[1:29 31:72]
for ii=1:length(conbinshi)
    avls.conbins{conbinshi(ii)}= [2400 2800]
end

for ii=1:length(conbinsmid)
    avls.conbins{conbinsmid(ii)}= [2270 2800]
end

for ii=1:length(conbinsnorm)
    avls.conbins{conbinsnorm(ii)}=[2100 2800]
end

    
    
    
    
    