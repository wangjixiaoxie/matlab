%master script, load and combine various fvs, from postscreen./50microstim,
%and 10 microstim

clear avls
avls.baspath='/oriole1/pk72r70'

pvl{1}='/templtest1';
cvl{1}='batch06.keep'
pvl{2}='/templtest2'
cvl{2}='batch.keep.rand'
pvl{3}='/templtest3'
cvl{3}='batch.keep.rand'
pvl{4}='/templtest4'
cvl{4}='batch.keep.rand'
pvl{5}='/templtest5'
cvl{5}='batch.keep.rand'
pvl{6}='/templtest8'
cvl{6}='batch'
pvl{7}='/preday'
cvl{7}='batch.keep.rand'


%start of WN feedback

%8-11-07 to 8-15-07
pvl{8}='/wnon'
cvl{8}='batchcomb.keep.catch'
%recovery?
%WN off  8-16 to 8-20???
%wnoff
pvl{9}='/wnoff/19files'
cvl{9}='batch.rand.keep.rand'



% second shift???
pvl{10}='/wnonsyn2'
cvl{10}='batch03.keep.catch'

pvl{11}='/wnonsyn2'
cvl{11}='batchcomb.keep.rand'

avls.pvls=pvl;
avls.cvl=cvl;


avls.bastms{1}={'2007-8-06' '2007-8-10'}
avls.wntms{1}={'2007-8-11' '2007-8-15'},
avls.rectms{1}={'2007-8-16' '2007-8-19'}
avls.offset_dys(1)=[4];
avls.bastms{2}={'2007-10-04' '2007-10-06'}
avls.wntms{2}={'2007-10-07' '2007-10-12'},
avls.rectms{2}={'2007-10-12' '2007-10-15'}
avls.offset_dys(2)=NaN;
avls.LEARNANALIND=[1 ]
avls.stabwin={'2007-08-06' '2007-08-10'}

avls.DOCONTPTANAL=0;
avls.DOCONTSYNANAL=0;

avls.pvls=pvl;
avls.cvl=cvl;
avls.ANALPTNT=1;
avls.SEQNT=1;
avls.DOCONT_ANAL=0;
for ii=1:length(pvl)
avls.SEC(ii)=0;
end
avls.fbins{1}=[3000 3500];
avls.fbins{2}=[2000 3000];
avls.fbins{3}=[2000 3000];
avls.fbins{4}=[2000 3000];


avls.NT{1}='e';
avls.NT{2}='a';
avls.NT{3}='i'
avls.NT{4}='c'
for kk=1:4
    avls.PRENT{kk}='d'
end
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025
avls.tbinshft(3)=.025;
avls.tbinshft(4)=.025;

avls.NFFT(1)=256;
avls.NFFT(2)=256;
avls.NFFT(3)=256;
avls.NFFT(4)=256;
for ii=1:length(avls.pvls)
    avls.EV4(ii)=0;
end



avls.MKFV=[1:11]
avls.PTIND=[]    
avls.SEQIND=1:11
% avls.bastms{1}={'2007-5-19' '2007-5-21'}
% avls.wntms{1}={'2007-5-22' '2007-5-30'},
% avls.rectms{1}={'2007-5-31' '2007-6-04'}
% avls.NTIND{1}=1;
% avls.NTIND{2}=2;
avls.ANALPTNT=3;
avls.SEQTRGNT{1}=2;
avls.ALLSEQNT=[1 2 3 4]

