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


avls.NT{1}='a';
avls.NT{2}='b';

avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025;
avls.fbins{1}=[3000 3500];
avls.fbins{2}=[2200 2700];
avls.NFFT(1)=256;
avls.NFFT(2)=512;
avls.MKFV=[10]
   