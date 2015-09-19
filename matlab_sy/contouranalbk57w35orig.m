
avls.pvls{1}='/oriole6/bk57w35/stimoff515/'
avls.cvl{1}='batch16.catch.keep.rand'
avls.NT{1}='a'

avls.pvls{2}='/oriole6/bk57w35/520_100ma_40ms_bilat2/'
avls.cvl{2}='batch'
avls.NT{2}='a'

avls.pvls{3}='/oriole6/bk57w35/521_100ma_40ms_bilat2'
avls.cvl{3}='batch'
avls.NT{3}='a'

avls.pvls{4}='/oriole6/bk57w35/522_100ma_40ms_400Hz/'
avls.cvl{4}='batch'
avls.NT{4}='a'

avls.pvls{5}='/oriole6/bk57w35/522_100ma_40ms_400Hz_2/'
avls.cvl{5}='batch'
avls.NT{5}='a'

avls.pvls{6}='/oriole6/bk57w35/523_100ma_400Hz_40ms/'
avls.cvl{6}='batch'
avls.NT{6}='a'

avls.pvls{7}='/oriole6/bk57w35/524_100mz_100Hz_40ms/'
avls.cvl{7}='batch'
avls.NT{7}='a'

avls.pvls{8}='/oriole6/bk57w35/524_100ma_400Hz_40ms/'
avls.cvl{8}='batch'
avls.NT{8}='a'

avls.pvls{9}='/oriole6/bk57w35/526_100mz_100Hz_40ms/'
avls.cvl{9}='batch26'
avls.NT{9}='a'

avls.pvls{10}='/oriole6/bk57w35/526_100ma_400Hz_40ms/'
avls.cvl{10}='batch'
avls.NT{10}='a'


avls.pvls{11}='/oriole6/bk57w35/526_100ma_40Hz_40ms/'
avls.cvl{11}='batch'
avls.NT{11}='a'

avls.pvls{12}='/oriole6/bk57w35/528_100ma_100Hz_40ms/'
avls.cvl{12}='batch28'
avls.NT{12}='a'

avls.pvls{13}='/oriole6/bk57w35/528_100ma_400Hz_40ms/'
avls.cvl{13}='batch28'
avls.NT{13}='a'

avls.pvls{14}='/oriole6/bk57w35/529_100Hz_100ma_newconfig/'
avls.cvl{14}='batch'
avls.NT{14}='a'

avls.pvls{15}='/oriole6/bk57w35/530_newconfig_100Hz_CORRECT/'
avls.cvl{15}='batch'
avls.NT{15}='a'

avls.pvls{16}='/oriole6/bk57w35/601_400Hz_60ms_oldconfig/'
avls.cvl{16}='batch'
avls.NT{16}='a'

avls.pvls{17}='/oriole6/bk57w35/602_100Hz_60ms_origconfig/'
avls.cvl{17}='batch'
avls.NT{17}='a'

avls.pvls{18}='/oriole6/bk57w35/602_400Hz_60ms_orig_2/'
avls.cvl{18}='batch'
avls.NT{18}='a'

avls.pvls{19}='/oriole6/bk57w35/602_40Hz_60ms_origconfig/'
avls.cvl{19}='batch02'
avls.NT{19}='a'

avls.pvls{20}='/oriole6/bk57w35/603_100Hz_60ms_origconfig/'
avls.cvl{20}='batch'
avls.NT{20}='a'

avls.pvls{21}='/oriole6/bk57w35/603_100Hz_origconfig_60ms/'
avls.cvl{21}='batch'
avls.NT{21}='a'

avls.pvls{22}='/oriole6/bk57w35/603_400Hz_60ms_origconfig/'
avls.cvl{22}='batch'
avls.NT{22}='a'

avls.pvls{23}='/oriole6/bk57w35/603_40Hz_60ms_origconfig/'
avls.cvl{23}='batch'
avls.NT{23}='a'


ptst.muinds=[]
% clear fvpt valsa
ptst.pt_NFFT=512;
ptst.con_NFFT=4096;
ptst.fbins=[6000 8000]
ptst.conbins=[2000 2500]
ptst.pt_tbinshft=.08;
ptst.con_tbinshft=.05;
ptst.avls=avls;

avls.SOUNDSTR='und'
avls.STIMSTR='im'

avls.STIMBNDS=[200 0]
avls.HST_EDGES=[6300:50:8000]




