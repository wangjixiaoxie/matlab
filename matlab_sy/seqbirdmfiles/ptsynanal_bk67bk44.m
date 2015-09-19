%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls
avls.baspath='/oriole7/dir2/bk67bk44/'


avls.pvls{1}='0915_wnonUP_FF/'
avls.cvl{1}='batch.catch'
avls.WINDOWSIZE=16;
avls.pvls{2}='0916pm_wnoff/'
avls.cvl{2}='batch.catch'


avls.pvls{3}='wnon_913am/13files/'
avls.cvl{3}='batch.catch.keep'

avls.pvls{4}='wnon_913am/14files/'
avls.cvl{4}='batch.catch.keep'

avls.pvls{5}='0924_wnoff/'
avls.cvl{5}='batch.rand'

avls.pvls{6}='postlesionsyn_925/'
avls.cvl{6}='batch.catch'

avls.pvls{7}='postlesion/'
avls.cvl{7}='batch21.rand'


avls.pvls{8}='postlesionwnon/'
avls.cvl{8}='batch.catch'


avls.pvls{9}='0927_wnoff/'
avls.cvl{9}='batch.catch'

avls.pvls{10}='0928_wnonpt/'
avls.cvl{10}='batch.catch'

avls.pvls{11}='wnon_913am/12files'
avls.cvl{11}='batch12.catch.keep'

avls.pvls{12}='synatemptest_910'
avls.cvl{12}='batch2.rand.keep'

avls.pvls{13}='wnoff_tmptest_914pm'
avls.cvl{13}='batch.catch'

avls.NT{1}='a';
avls.NT{2}='b';
for ii=1:2
    avls.PRENT{ii}='';
end
for ii=1:length(avls.pvls)
avls.EV4(ii)=0;
end
%sequence pre lesion
avls.bastms{1}={'2010-9-11' '2010-9-12'}
avls.wntms{1}={'2010-9-13' '2010-9-14'}
avls.rectms{1}={'2010-9-15' '2010-9-15'}
avls.type{1}='SQ'
avls.offset_dys(1)=1;
avls.SPLITDAY(2)=1;
%pitch prelesion

avls.bastms{2}={'2010-9-15' '2010-9-15' }
avls.wntms{2}={'2010-9-16' '2010-9-16'}
avls.rectms{2}={'2010-9-17' '2010-9-17'}
avls.type{2}='PT'
avls.drxn{2}='up'
avls.offset_dys(2)=0;
avls.SPLITDAY(2)=0;
%pitch post lesion(2)
avls.bastms{3}={'2010-9-28' '2010-9-28' }
avls.wntms{3}={'2010-9-29' '2010-9-30'}
avls.rectms{3}={'2010-10-01' '2010-10-01'}
avls.type{3}='PT'
avls.drxn{3}='up'
avls.offset_dys(3)=1;
avls.SPLITDAY(3)=0;


%pitch post lesion(1)
avls.bastms{4}={'2010-9-21' '2010-9-21'}
avls.wntms{4}={'2010-9-22' '2010-9-23'}
avls.rectms{4}={'2010-9-24' '2010-9-26'}
avls.type{4}='PT'
avls.drxn{4}='up'
avls.offset_dys(4)=1;
avls.SPLITDAY(4)=1;
%sequence post lesion(1)
avls.bastms{5}={'2010-9-23' '2010-9-24' }
avls.wntms{5}={'2010-9-25' '2010-9-26'}
avls.rectms{5}={'2010-9-27' '2010-9-27' }
avls.type{5}='SQ'
avls.offset_dys(5)=1;
avls.SPLITDAY(5)=0;

avls.seq_days{1}={'2010-9-11' '2010-9-16'}
avls.seq_start(1)=2;
avls.seq_days{2}={'2010-9-23' '2010-9-28'}
avls.seq_start(2)=2;
avls.pt_days{1}={'2010-9-13' '2010-9-16'}
avls.pt_start(1)=3;
avls.pt_days{2}={'2010-9-28' '2010-10-01'}
avls.pt_start(2)=1;

avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025;
avls.fbins{1}=[3000 3500];
avls.fbins{2}=[2200 2700];
avls.NFFT(1)=256;
avls.NFFT(2)=512;
avls.MKFV=[1:13]
avls.PTIND=[]
avls.DOCONTSYNANAL=0;
avls.DOCONTPTANAL=0;
for ii=1:5
    avls.SEQTRGNT{ii}=1;
end
avls.ALLSEQNT=[1:2]
% avls.SPLITDAYS=1;
avls.ANALPTNT=2;
avls.SEQIND=[1:13];
avls.LESIONDATE={'2010-9-18'}
avls.SPLITDAY=1;



   