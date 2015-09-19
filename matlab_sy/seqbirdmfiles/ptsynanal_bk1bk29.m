%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2

%why two pitch shifts??

clear fv fbins bt smoothamp baspath conbins avls
avls.baspath='/oriole7/dir1/bk1bk29/'
avls.pvls{1}='syntest_B/'
avls.cvl{1}='batchcomb'
avls.WINDOWSIZE=16;
avls.pvls{2}='0913_wnon/'
avls.cvl{2}='batch.catch.keep'

avls.pvls{3}='0915_wnoff/'
avls.cvl{3}='batch.catch'

avls.pvls{4}='0922_pitch_wnon/'
avls.cvl{4}='batch.catch'

avls.pvls{5}='0916_wnon/'
avls.cvl{5}='batch.catch'


avls.pvls{6}='0916pm_wnoff/'
avls.cvl{6}='batch.catch'

avls.pvls{7}='postlesion_pt_test/'
avls.cvl{7}='batch.catch'

avls.pvls{8}='0922_pitch_wnon/'
avls.cvl{8}='batch.catch'


avls.pvls{8}='0924_wnoff/'
avls.cvl{8}='batch.catch'

avls.pvls{9}='0927_wnoff/'
avls.cvl{9}='batch.catch'

avls.pvls{10}='prept_test_0928/'
avls.cvl{10}='batch.catch'

avls.pvls{11}='postlesion/'
avls.cvl{11}='batch20.catch'

avls.pvls{12}='syntest_B/'
avls.cvl{12}='batch11.catch'



avls.NT{1}='a';
avls.NT{2}='b';
avls.NT{3}='c'
for ii=1:3
    avls.PRENT{ii}='';
end
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025
avls.tbinshft(3)=.025;

for ii=1:length(avls.pvls)
avls.EV4(ii)=0;
end



avls.fbins{1}=[3000 3500];
avls.fbins{2}=[2000 3000];
avls.fbins{3}=[3000 3500];


%sequence pre lesion
avls.bastms{1}={'2010-9-11' '2010-9-12'}
avls.wntms{1}={'2010-9-13' '2010-9-14'}
avls.rectms{1}={'2010-9-15' '2010-9-15'}
avls.type{1}='SQ'
avls.offset_dys(1)=1;
%pitch prelesion

avls.bastms{2}={'2010-9-15' '2010-9-15' }
avls.wntms{2}={'2010-9-16' '2010-9-16'}
avls.rectms{2}={'2010-9-17' '2010-9-17'}
avls.type{2}='PT'
avls.drxn{2}='up'
avls.offset_dys(2)=0;
%pitch post lesion(2)
avls.bastms{3}={'2010-9-28' '2010-9-28' }
avls.wntms{3}={'2010-9-29' '2010-9-30'}
avls.rectms{3}={'2010-10-01' '2010-10-01'}
avls.type{3}='PT'
avls.drxn{3}='up'
avls.offset_dys(3)=1;
%pitch post lesion(1)
avls.bastms{4}={'2010-9-20' '2010-9-21'}
avls.wntms{4}={'2010-9-22' '2010-9-23'}
avls.rectms{4}={'2010-9-24' '2010-9-26'}
avls.type{4}='PT'
avls.drxn{4}='up'
avls.offset_dys(4)=1;
% sequence post lesion(1)
avls.bastms{5}={'2010-9-24' '2010-9-24' }
avls.wntms{5}={'2010-9-25' '2010-9-26'}
avls.rectms{5}={'2010-9-27' '2010-9-27' }
avls.type{5}='SQ'
avls.offset_dys(5)=1;

avls.seq_days{1}={'2010-9-11' '2010-9-16'}
avls.seq_start(1)=2
avls.seq_days{2}={'2010-9-23' '2010-9-28'}
avls.seq_start(2)=2
avls.pt_days{1}={'2010-9-13' '2010-9-17'}
avls.pt_start(1)=3
avls.pt_days{2}={'2010-9-20' '2010-9-24'}
avls.pt_start(2)=2;

avls.ANALPTNT=[1 3];
avls.LESIONDATE={'2010-9-18'}


avls.NFFT(1)=256;
avls.NFFT(2)=512;
avls.NFFT(3)=256;
avls.MKFV=[1:12]
avls.PTIND=[1:12]
avls.SEQIND=[1:12]
for ii=1:5
    avls.SEQTRGNT{ii}=2
end
avls.ALLSEQNT=[1:2]
avls.DOCONTSYNANAL=0;
avls.DOCONTPTANAL=0;

avls.SPLITDAY=1;
