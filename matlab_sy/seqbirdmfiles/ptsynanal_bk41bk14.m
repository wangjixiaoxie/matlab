;%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls
avls.baspath='/oriole2/bk41bk14/'

%5-02 and 5-03
avls.pvls{1}='temptest2/'
avls.cvl{1}='batchcomb'
avls.WINDOWSIZE=16;

avls.pvls{2}='temptest3/'
avls.cvl{2}='batch04.catch.keep'

% avls.pvls{3}='syntest3/'
% avls.cvl{3}='batch.catch'

avls.pvls{3}='wnonsyn1/'
avls.cvl{3}='batchcomb'

avls.pvls{4}='wnoff511_840/'
avls.cvl{4}='batchcomb'

avls.pvls{5}='518comb/'
avls.cvl{5}='batch18.catch'

avls.pvls{6}='astim519/'
avls.cvl{6}='batchcomb'




avls.pvls{7}='70109data/'
avls.cvl{7}='batch01.catch'

avls.pvls{8}='wnonsyn3/'
avls.cvl{8}='batch03.catch'

avls.pvls{9}='070409/'
avls.cvl{9}='batch04.catch'

avls.pvls{10}='wnonsyn3/'
avls.cvl{10}='batch05.catch'

avls.pvls{11}='728/'
avls.cvl{11}='batch28.catch'

avls.pvls{12}='wnapitchpstlesion/'
avls.cvl{12}='batch.catch'

avls.pvls{13}='screen/'
avls.cvl{13}='batch30'
avls.pvls{14}='temptest2/'
avls.cvl{14}='batch01.keep.rand'
avls.pvls{15}='wnoff511_840/11files'
avls.cvl{15}='batch.catch'
avls.pvls{16}='wnonsyn1'
avls.cvl{16}='batch09.catch'



% avls.pvls{7}='wnoff511_840/'
% avls.cvl{7}='batchcomb'



%prelesion sequence shift
avls.bastms{1}={'2009-05-04' '2009-05-04'}
avls.wntms{1}={'2009-05-05' '2009-05-10'},
avls.rectms{1}={'2009-05-11' '2009-05-13'}
avls.type{1}='SQ'
avls.offset_dys(1)=3;


%prelesion pitch shift
avls.bastms{2}={'2009-05-18' '2009-05-18'}
avls.wntms{2}={'2009-05-19' '2009-05-21'}
avls.rectms{2}={}
avls.type{2}='PT'
avls.drxn{2}='up'
avls.offset_dys(2)=2;

%postlesion sequence shift
avls.bastms{3}={'2009-07-01' '2009-07-01'}
avls.wntms{3}={'2009-07-03' '2009-07-05'}
avls.rectms{3}={'2009-07-06' '2009-07-07'}
avls.type{3}='SQ'
avls.offset_dys(3)=1


%postlesion pitch
avls.bastms{4}={'2009-07-28' '2009-07-28'}
avls.wntms{4}={'2009-08-05' '2009-08-07'}
avls.rectms{4}=[];
avls.type{4}='PT'
avls.drxn{4}='up'
avls.offset_dys(4)=2


avls.stabwin={'2009-05-01' '2009-05-04'}
avls.LEARNANALIND=1


%postlesion pitch shift  %no matched data!!
% avls.bastms{4}={'2011-01-10' '2011-01-11'}
% avls.wntms{4}={'2011-01-12' '2011-01-14'}
% avls.rectms{4}={'2011-01-15' '2011-01-16'}


avls.DOCONTPTANAL=0;
avls.DOCONTSYNANAL=0;


avls.NT{1}='a';
avls.NT{2}='b';
for ii=1:2
    avls.PRENT{ii}='';
end
for ii=1:length(avls.pvls)
    avls.EV4(ii)=0;
end
% avls.NT{3}='c'
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025
% avls.tbinshft(3)=.025;


avls.fbins{1}=[3000 3500];
avls.fbins{2}=[2000 3000];
% avls.fbins{3}=[3000 3500];


avls.ANALPTNT=1;
avls.SEQTRGNT{1}=1;
avls.ALLSEQNT=[1:2]
avls.SPLITDAYS=0;


avls.NFFT(1)=256;
avls.NFFT(2)=512;
% avls.NFFT(3)=256;
avls.MKFV=[1:16]
avls.PTIND=[16]
avls.SEQIND=1:16
avls.LESIONDATE={'2009-5-30'}
avls.SPLITDAY=1
   