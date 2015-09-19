%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls
avls.baspath='/oriole7/dir2/r29g89/'

avls.pvls{1}='1109_tmptest2/'
avls.cvl{1}='batch10.catch'


avls.pvls{2}='1111_wnfixed/'
avls.cvl{2}='batch12.catch.keep'
avls.WINDOWSIZE=16;
avls.pvls{3}='1111_wnfixed/13files/'
avls.cvl{3}='batch.catch.keep'

avls.pvls{4}='1111_wnfixed/14files/'
avls.cvl{4}='batch14.catch'

avls.pvls{5}='1111_wnfixed/15files/'
avls.cvl{5}='batch.catch'
avls.pvls{6}='1111_wnfixed/16files/'
avls.cvl{6}='batch16.catch'

% avls.pvls{5}='1130_pttest/'
% avls.cvl{5}='batch01.catch'

avls.pvls{7}='1208_prelesionwnoff/'
avls.cvl{7}='batch11.catch'

avls.pvls{8}='1208_prelesionwnoff/'
avls.cvl{8}='batch12.catch'

avls.pvls{9}='1212_wnon/12files'
avls.cvl{9}='batch.catch'

avls.pvls{10}='1212_wnon/13files'
avls.cvl{10}='batch.catch'

avls.pvls{11}='1212_wnon/14files/'
avls.cvl{11}='batch.catch'


%pitch shifts
avls.pvls{12}='1130_pttest/'
avls.cvl{12}='batch01.catch'

avls.pvls{13}='1202wnonprelesion/'
avls.cvl{13}='batch.catch'

avls.pvls{14}='1218postlesion_2/'
avls.cvl{14}='batch.catch'

avls.pvls{15}='postlesion/0104wnon'
avls.cvl{15}='batch.catch'

avls.pvls{16}='postlesion/0104wnon/06files'
avls.cvl{16}='batch06.catch'

avls.pvls{17}='postlesion/0104wnon/07files'
avls.cvl{17}='batch.rand'

avls.pvls{18}='postlesion/0104wnon'
avls.cvl{18}='batch10.catch'

avls.pvls{19}='postlesion/0104wnon'
avls.cvl{19}='batch10.catch'

avls.pvls{20}='postlesion/0111_synconfigtest'
avls.cvl{20}='batch11.catch'

avls.pvls{21}='postlesion/011210_wnon'
avls.cvl{21}='batchcomb'
avls.pvls{22}='postlesion/011510_wnoff'
avls.cvl{22}='batchcomb'

avls.pvls{23}='postlesion/0104wnon/09files'
avls.cvl{23}='batch09.rand'

for ii=1:length(avls.pvls)
    avls.EV4(ii)=0
end

%prelesion sequence shift
avls.bastms{1}={'2010-11-10' '2010-11-10'}
avls.wntms{1}={'2010-11-12' '2010-11-13'},
avls.rectms{1}={'2010-11-14' '2010-11-15'}
avls.type{1}='SQ'
avls.offset_dys(1)=1;

%prelesion pitch shift
avls.bastms{2}={'2010-12-11' '2010-12-11'}
avls.wntms{2}={'2010-12-12' '2010-12-14'}
avls.rectms{2}={}
avls.type{2}='PT'
avls.drxn{2}='dn'
avls.offset_dys(2)=1;

%postlesion pitch shift
avls.bastms{3}={'2011-01-03' '2011-01-03'}
avls.wntms{3}={'2011-01-04' '2011-01-06'}
avls.rectms{3}={'2011-01-09' '2011-01-09'}
avls.type{3}='PT'
avls.drxn{3}='dn'
avls.offset_dys(3)=1;

%postlesion sequence shift
avls.bastms{4}={'2011-01-10' '2011-01-11'}
avls.wntms{4}={'2011-01-12' '2011-01-14'}
avls.rectms{4}={'2011-01-15' '2011-01-16'}
avls.offset_dys(4)=1;
avls.type{4}='SQ'
avls.MKFV=[1:23]
avls.SEQIND=[1:23]
avls.PTIND=[]

avls.seq_days{1}={'2010-11-10' '2010-11-15'}
avls.seq_start(1)=2;
avls.seq_days{2}={'2011-1-10' '2011-1-16'}
avls.seq_start(2)=2;
avls.pt_days{1}={'2010-12-11' '2010-12-14'}
avls.pt_start(1)=1;
avls.pt_days{2}={'2011-1-03' '2011-1-09'}
avls.pt_start(2)=1;


for ii=1:length(avls.pvls)

   avls.SEC(ii)=0; 
end
avls.NT{1}='b';
avls.NT{2}='a';
avls.NT{3}='c'
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.028
avls.tbinshft(3)=.025;

for ii=1:3
    avls.PRENT{ii}=''
end
for ii=1:22
    avls.EV4(ii)=0
end

avls.DOCONT_ANAL=1;
avls.fbins{1}=[3000 3500];
avls.fbins{2}=[3000 3500];
avls.fbins{3}=[3000 3500];


avls.ANALPTNT=2;
for ii=1:4
avls.SEQTRGNT{ii}=2;
end
avls.ALLSEQNT=[1 2 3]
avls.SPLITDAYS=0;


avls.NFFT(1)=256;
avls.NFFT(2)=256;
avls.NFFT(3)=256;


avls.DOCONTPTANAL=0;
avls.DOCONTSYNANAL=0;
avls.SPLITDAY=1;
avls.LESIONDATE={'2010-12-20'}