%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls


avls.baspath='/oriole5/r3g7/'
avls.datdir='datasum'

wnin(1).tmon{1}={'2009-8-22 7' '2009-8-24 7'}
wnin(1).tmon{2}={'2009-9-02 7' '2009-9-04 7'}
wnin(1).tmon{3}={'2009-9-26 7'}
wnin(1).tmon{4}={'2009-10-31 7'}
wnin(1).tmoff{1}={'2009-8-24 7' '2009-8-29 7'}
wnin(1).tmoff{2}={'2009-9-04 7' '2009-9-09 7'} 
wnin(1).tmoff{3}={'2009-10-13 7'}
wnin(1).tmoff{4}={'2009-11-10 7'}


    wnrevin(1).tmon{1}={'2009-8-29 7'}
wnrevin(1).tmoff{1}={'2009-9-02 7'}
wnrevin(1).tmon{2}={'2009-9-09 7'}
wnrevin(1).tmoff{2}={'2009-9-12 7'}
wnrevin(1).tmon{2}={'2009-9-09 7'}
wnrevin(1).tmoff{2}={'2009-9-12 7'}
wnrevin(1).tmon{3}={'2009-10-13 7'}
wnrevin(1).tmoff{3}={'2009-10-16 7'}
wnrevin(1).tmon{4}={'2009-11-11 7'}
wnrevin(1).tmoff{4}={'2009-11-17 7'}


avls.wnin=wnin
avls.wnrevin=wnrevin



%leadin
%temptest_aug21
%temptest

avls.pvls{1}='0427_acsf/'
avls.cvl{1}='batch'
avls.NT{1}='c'

avls.pvls{2}='0428_ap5am'
avls.cvl{2}='batch'
avls.NT{2}='c'
avls.pvls{3}='0424_acsf_wnoff'
avls.cvl{3}='batch.keep'
avls.NT{3}='c'

avls.pvls{4}='0427_5mMapv'
avls.cvl{4}='batch'
avls.NT{4}='c'

avls.pvls{5}='0324_probesin'
avls.cvl{5}='batch.keep.rand'
avls.NT{5}='c'

avls.pvls{6}='0326_5mMapv'
avls.cvl{6}='batch.keep'
avls.NT{6}='c'

avls.NOSTIM(1)=0;
for ii=1:8
    
    avls.del{ii}=.065
end
for ii=9
    avls.del{ii}=.095;
end
%do pitch analysis on these batch inds.

avls.mkfv=[5:6]

    avls.con_tbinshft=0;
    %shift for pitch analysis
%     avls.pt_tbinshft=.072;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  %shift for contours
   for ii=1:length(avls.pvls)
    avls.conbins{ii}=[5600 7200]
   end
    avls.fbins=[3000 3400]
%     avls.conbins=[2400 2800];
    
    
    avls.basruns=[]
    avls.contanal=[5:6]
avls.contfv=[5:6]
avls.ptfv=[5:6]
avls.analfv=[5:6]
avls.catchstimfv=[2]
avls.contnum=[1 ];
    avls.con_tbinshft=-.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.074;
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
    avls.STAN_RUNS=[1:7 13:15 17:19 23:25 28:33 35:36 39 41:43 45 52 54:55 57:63 65:68 70 75 77:87]
    avls.REMOVEOUTLIERS=1;
    avls.NUMSTD=2;
    conbinsmid=[30 70]	
    conbinshi=[73:87]
    conbinsnorm=[1:29 31:47 52:69 71:72 88:92]
%     conbinspec=[48 49 50 51]
% for ii=1:length(conbinshi)
%     avls.conbins{conbinshi(ii)}= [2400 2800]
% end
% 
% for ii=1:length(conbinsmid)
%     avls.conbins{conbinsmid(ii)}= [2270 2800]
% end
% 
% for ii=1:length(conbinsnorm)
%     avls.conbins{conbinsnorm(ii)}=[2100 2800]
% end
% 
% for ii=1:length(conbinspec)
%     avls.conbins{conbinspec(ii)}=[2150 2450]
% end
%     
%     
%     
%     