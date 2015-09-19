%this script is written to generate pitch contours for bk61w42.
clear fv fbins bt smoothamp baspath conbins

% avls.pvls{1}='/oriole6/bk57w35/test9/'
% avls.cvl{1}='batch.notcatch'
% avls.NT{1}='a'

% tvls.pvls{1}='/oriole6/bk57w35/617_125ma_100Hz_40ms/'

baspath='/oriole/bk48w74/'
% 
% 
% %short delay
% out1_2350ma_out2_1440ma_90msdel_20ms_200Hz


% [01;34mstim918_100msdel_70ma_20msburst[0m
% [01;34mstimon718_100msdel_20msburst_90ma[0m
% [01;34mstimon718_90msdel_10msburst_70maout1_23_out2_14[0m
% [01;34mstimon718_90msdel_10msburst_90maout1_23_out2_14[0m
% 
% 
% stim719_110msdel_10ms_80ma
% [01;34mstim719_3_100msdel_80ma_10mspulse[0m
% 
% %720
% [0m[01;34mstim720_3mspulse_200Hz_80ma_105msdel[0m
% [01;34mstimon720_5mspulse_80ma_200Hz[0m
% [01;34mstimon720_5mspulse_80ma_200Hz_105msdel[0m
% 
% [01;34mstim722_5ms_85mdel_70ma_200Hz[0m
% 
% %725
% stim725_85msdel_80ma_5ms_200Hz
% 
% %initial baseline 714_files

tvls.pvls{1}='stim14_bilat_70ma_65msdel_200Hz/'
tvls.cvl{1}='batch'
tvls.NT{1}='a'

tvls.pvls{2}='stim14_bilat_70ma_65msdel_200Hz_run2/'
tvls.cvl{2}='batch'
tvls.NT{2}='a'

tvls.pvls{3}='stim23_bilat_20ma_65msdel_200Hz_aug20/'
tvls.cvl{3}='batch'
tvls.NT{3}='a'

tvls.pvls{4}='stim14_bilat_90ma_65msdel_200Hz/'
tvls.cvl{4}='batch'
tvls.NT{4}='a'


tvls.pvls{5}='stimaug23/'
tvls.cvl{5}='batch'
tvls.NT{5}='a'


% tvls.pvls{2}='714/out1_23_40ma_out2_1450ma_200Hz_70msdel/'
% tvls.cvl{2}='batch'
% tvls.NT{2}='a'
% 
% tvls.pvls{3}='714/out1_23_40ma_out2_1450ma_200Hz_70msdel/'
% tvls.cvl{3}='batch'
% tvls.NT{3}='a'
% 
% 
% %716files -- first day of pitchshift
% tvls.pvls{4}='716/stimon716_35maboth'
% tvls.cvl{4}='batch'
% tvls.NT{4}='a'
% 
% tvls.pvls{5}='716/stimon716_40maout1_35maout2'
% tvls.cvl{5}='batch'
% tvls.NT{5}='a'
% 
% tvls.pvls{6}='716/stimon716_40maout1_35maout2'
% tvls.cvl{6}='batch'
% tvls.NT{6}='a'
% 
% %717files  -- first day of pitchshift
% tvls.pvls{7}='717/stim71713out2_14out1_50ma'
% tvls.cvl{7}='batch'
% tvls.NT{7}='a'
% 
% tvls.pvls{8}='717/stimon717_13out2_24_out1'
% tvls.cvl{8}='batch'
% tvls.NT{8}='a'
% 
% tvls.pvls{9}='717/stimon717_23out2_14out1_40ma_60ms_75del_200Hz'
% tvls.cvl{9}='batch'
% tvls.NT{9}='a'
% 
% %718files 
% 
% tvls.pvls{10}='718/stim718_400Hz_40ms_75msdel'
% tvls.cvl{10}='batch'
% tvls.NT{10}='a'
% 
% 
% tvls.pvls{11}='718/stimon718_40ms_75del_out2_14_out1_23'
% tvls.cvl{11}='batch'
% tvls.NT{11}='a'
% 
% %719 files
% 
% tvls.pvls{12}='719/stim719_overnight_75msdel_60mspulse_200Hz'
% tvls.cvl{12}='batch'
% tvls.NT{12}='a'
% 
% tvls.pvls{13}='719/stimprobe719'
% tvls.cvl{13}='batch'
% tvls.NT{13}='a'
% 
% tvls.pvls{14}='/719/stim719-2/'
% tvls.cvl{14}='batch'
% tvls.NT{14}='a'
% 
% 
% %720 files
% 
% tvls.pvls{15}='720/stimon720_60ms_200Hz'
% tvls.cvl{15}='batch'
% tvls.NT{15}='a'
% 
% %721files
% tvls.pvls{16}='721/stim721_STIMONCOR_75del_200Hz_40ma_60ms'
% tvls.cvl{16}='batch'
% tvls.NT{16}='a'
% 
% %722files
% tvls.pvls{17}='722/stim722_40ms_40ma_75msdel'
% tvls.cvl{17}='batch'
% tvls.NT{17}='a'
% 
% tvls.pvls{18}='722/stimon_standard_60ms_200Hz_75del_722'
% tvls.cvl{18}='batch'
% tvls.NT{18}='a'
% 
% %723files
% 
% tvls.pvls{19}='723/stim723_75msdel_40ms_200Hz_40ma'
% tvls.cvl{19}='batch'
% tvls.NT{19}='a'
% 
% %724files???what happened 7/24
% tvls.pvls{20}='724/bilatstim724_45ma_40ms_75msdel'
% tvls.cvl{20}='batch'
% tvls.NT{20}='a'
% 
% %725 files
% tvls.pvls{21}='725/stim725_50ma_40ms_200Hz'
% tvls.cvl{21}='batch'
% tvls.NT{21}='a'
% 
% %726files
% tvls.pvls{22}='726/stim726_50ma_40ms_75del'
% tvls.cvl{22}='batch'
% tvls.NT{22}='a'
% 
% %727files
% tvls.pvls{23}='727/stim727_standard_50ma'
% tvls.cvl{23}='batch'
% tvls.NT{23}='a'
% 
% tvls.pvls{24}='728/stim728_standard'
% tvls.cvl{24}='batch'
% tvls.NT{24}='a'
% 
% 
% tvls.pvls{25}='730/stim730'
% tvls.cvl{25}='batch'
% tvls.NT{25}='a'
% 
% tvls.pvls{26}='731/stim731'
% tvls.cvl{26}='batch.rand'
% tvls.NT{26}='a'
% 
% tvls.pvls{27}='801/stim801'
% tvls.cvl{27}='batch'
% tvls.NT{27}='a'
% 
% tvls.pvls{28}='802/stim802'
% tvls.cvl{28}='batch'
% tvls.NT{28}='a'
% 
% % tvls.pvls{25}='wn729rev'
% % tvls.cvl{25}='batch'
% % tvls.NT{25}='a'
% % 
% % tvls.pvls{26}='stimstandard_test_729'
% % tvls.cvl{26}='batch'
% % tvls.NT{26}='a'
% % tvls.pvls{6}='/oriole6/bk59w37redo/stimon718_90msdel_10msburst_70maout1_23_out2_14'
% % % tvls.cvl{6}='batch'
% % % tvls.NT{6}='a'
% % % tvls.stim{6}=[0.075 0.105]
% % 
% % tvls.pvls{7}='/oriole6/bk59w37redo/wnrev719'
% % tvls.cvl{7}='batch'
% % tvls.NT{7}='a'
% % tvls.stim{7}=[0.06 0.08]
% % 
% % 
% % tvls.pvls{8}='/oriole6/bk59w37redo/stim718_400Hz_40ms_75msdel'
% % tvls.cvl{8}='batch'
% % tvls.NT{8}='a'
% % tvls.stim{8}=[0.07 0.10]
% % %original pitchshifted 400Hz runs
% % 
% % 
% % tvls.pvls{7}='/oriole6/bk59w37redo/stimprobe719'
% % tvls.cvl{7}='batch'
% % tvls.NT{7}='a'
% % 
% % tvls.pvls{8}='/oriole6/bk59w37redo/stim719-2'
% % tvls.cvl{8}='batch'
% % tvls.NT{8}='a'
% % 
% % tvls.pvls{11}='/oriole6/bk59w37redo/stim719_3_100msdel_80ma_10mspulse'
% % tvls.cvl{11}='batch'
% % tvls.NT{11}='a'
% % 
% % tvls.pvls{12}='/oriole6/bk59w37redo/stim719_110msdel_10ms_80ma'
% % tvls.cvl{12}='batch'
% % tvls.NT{12}='a'
% % 
% % tvls.pvls{13}='/oriole6/bk59w37redo/stim719_overnight_75msdel_60mspulse_200Hz'
% % tvls.cvl{13}='batcheve'
% % tvls.NT{13}='a'
% % 
% % tvls.pvls{14}='/oriole6/bk59w37redo/stim719_overnight_75msdel_60mspulse_200Hz'
% % tvls.cvl{14}='batchmorn'
% % tvls.NT{14}='a'
% % 
% % tvls.pvls{15}='/oriole6/bk59w37redo/stimon720_60ms_200Hz'
% % tvls.cvl{15}='batch'
% % tvls.NT{15}='a'
% % 
% % % tvls.pvls{15}='/oriole6/bk57w35/602_400Hz_60ms_origconfig'
% % % tvls.cvl{15}='batch'
% % % tvls.NT{15}='a'
% % 
% % tvls.pvls{16}='/oriole6/bk59w37redo/stimon720_5mspulse_80ma_200Hz_105msdel'
% % tvls.cvl{16}='batch'
% % tvls.NT{16}='a'
% % 
% % % %2nd round of pitchshifted runs
% % % tvls.pvls{17}='/oriole6/bk59w37redo/stim720_3mspulse_200Hz_80ma_105msdel'
% % % tvls.cvl{17}='batch'
% % % tvls.NT{17}='a'
% % 
% % % tvls.pvls{18}='/oriole6/bk59w37redo/stim721_STIMONCOR_75del_200Hz_40ma_60ms'
% % % tvls.cvl{18}='batch'
% % % tvls.NT{18}='a'
% % 
% % tvls.pvls{19}='/oriole6/bk59w37redo/stimon_standard_60ms_200Hz_75del_722'
% % tvls.cvl{19}='batch'
% % tvls.NT{19}='a'
% % 
% % % tvls.pvls{20}='/oriole6/bk59w37redo/stim722_5ms_85mdel_70ma_200Hz'
% % % tvls.cvl{20}='batch'
% % % tvls.NT{20}='a'
% % 
% % % tvls.pvls{21}='/oriole6/bk59w37redo/stim722_40ms_40ma_75msdel'
% % % tvls.cvl{21}='batch'
% % % tvls.NT{21}='a'
% % 
% % % tvls.pvls{22}='/oriole6/bk59w37redo/stim723_75msdel_40ms_200Hz_40ma'
% % % tvls.cvl{22}='batch'
% % % tvls.NT{22}='a'
% % 
% % 
% % tvls.pvls{23}='/oriole6/bk59w37redo/stim723_75msdel_40ms_200Hz_40ma'
% % tvls.cvl{23}='batch24'
% % tvls.NT{23}='a'
% % tvls.modoffset{23}=[];
% % 
% % tvls.pvls{24}='/oriole6/bk59w37redo/stim725_50ma_40ms_200Hz'
% % tvls.cvl{24}='batch'
% % tvls.NT{24}='a'
% % tvls.modoffset{24}=[];
% % 
% % tvls.pvls{25}='/oriole6/bk59w37redo/stim725_85msdel_80ma_5ms_200Hz'
% % tvls.cvl{25}='batch'
% % tvls.NT{25}='a'
% % tvls.modoffset{25}=[];
% % 
% % tvls.pvls{26}='/oriole6/bk59w37redo/stim726_50ma_40ms_75del'
% % tvls.cvl{26}='batch'
% % tvls.NT{26}='a'
% % tvls.modoffset{26}=[];
% % 
% % tvls.pvls{27}='/oriole6/bk59w37/bilat_34_400Hz_150ma_60ms_zerodel'
% % tvls.cvl{27}='batch'
% % tvls.NT{27}='a'
% % tvls.modoffset{27}=[];
% % 
% % tvls.pvls{28}='/oriole6/bk59w37/bilat_23_20ms_125ma_20msdel_400Hz'
% % tvls.cvl{28}='batch'
% % tvls.NT{28}='a'
% tvls.modoffset{28}=[]
% 
% tvls.pvls{29}='/oriole6/bk59w37/bilat_12_400Hz_125ma_20ms_20msdel'
% tvls.cvl{29}='batch'
% tvls.NT{29}='a'
% tvls.modoffset{29}=[]
% 
% tvls.pvls{30}='/oriole6/bk59w37/bilat_12_20ms_10msdel_400Hz_125ma'
% tvls.cvl{30}='batch'
% tvls.NT{30}='a'
% tvls.modoffset{30}=[]
% 
% tvls.pvls{31}='/oriole6/r34w20/wnon726'
% tvls.cvl{31}='batchlim'
% tvls.NT{31}='b'
% tvls.modoffset{31}=[]
% 
% tvls.pvls{32}='/oriole6/bk59w37/implantin'
% tvls.cvl{32}='batchpre'
% tvls.NT{32}='a'
% tvls.modoffset{32}=[]
% 
% 
% tvls.pvls{33}='/oriole6/bk59w37/705_bilat12_125ma_60ms_zerodel'
% tvls.cvl{33}='batch'
% tvls.NT{33}='a'
% tvls.modoffset{33}=[]
% 
% tvls.pvls{34}='/oriole6/bk59w37/705_bilat_23_125ma_60ms_zerodel'
% tvls.cvl{34}='batch'
% tvls.NT{34}='a'
% tvls.modoffset{34}=[]
% 
% 
% % tvls.pvls{19}='/oriole6/bk57w35/603_400Hz_60ms_origconfig'
% % tvls.cvl{19}='batch'
% % tvls.NT{19}='a'
% % 
% % tvls.pvls{20}='/oriole6/bk57w35/603_400Hz_60ms_origconfig'
% % tvls.cvl{20}='batch'
% % tvls.NT{20}='a'
% % 
% % 
% % tvls.pvls{21}='/oriole6/bk57w35/603_400Hz_60ms_origconfig'
% % tvls.cvl{21}='batch'
% % tvls.NT{21}='a'
% % 
% 
% 
% 


% muinds=[1 4 9:23]

% muinds=[1 4 9:24]
muinds=[5]
for ii=1:length(muinds)
    crind=muinds(ii)
    pathvl=tvls.pvls{crind}
    
    if(exist('baspath'))
        cmd=['cd ' baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end
    eval(cmd);
%    load sumdata.mat
    

    con_tbinshft=.05;
    pt_tbinshft=.079;
    pt_NFFT=512
    con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    
    
            fbins{crind}=[6000 8000]
            conbins{crind}=[2050 2700];
      
    
   
    NT=tvls.NT{crind};PRENT='';PSTNT='';
    bt{crind}=tvls.cvl{crind};
    fvpt=findwnote4(bt{crind},NT, PRENT, PSTNT, pt_tbinshft,fbins{crind}, pt_NFFT,1,'obs0');
    fv=findwnote4(bt{crind},NT,PRENT,PSTNT,con_tbinshft,fbins{crind},con_NFFT,1,'obs0');
    pitchdata{crind}=jc_pitchcontourFV(fv,1024,1020,1,conbins{crind}(1),conbins{crind}(2),[3 ],'obs0');
    contours=pitchdata{crind};
%     for ii=1:length(fv)
%        tmp=evsmooth(fv(ii).datt,32000)
%        smoothamp{crind}(:,ii)=resample(tmp,length(pitchdata{crind}(:,1)),length(tmp))
%     end
    tms=maketimevec2(bt{crind})
    cmd=['save ' bt{crind} '.mat contours tms fvpt'];
    eval(cmd);
end
 
initsamptime=512/32000;
initsamptimediff=1/8000;
pitchtms=initsamptime:initsamptimediff:(4096-512)/32000
pitchtms=pitchtms+con_tbinshft
figure
subplot(211)
%get the right chunk of data to make example plot
% 
% e[sm,sp,t,f]=evsmooth(fv{3}(58).datt,32000,10,512,0.8,2,100,10000);
% imagesc(t+.016,f,log(abs(sp)));syn;ylim([0,1e4]);
% hold on;
% plot(pitchtms,pitchdata{3}(:,58),'c')
% 
% 
% 
% subplot(212)
% 
% %generate 3 random pre
% inds{1}=[floor(rand(10,1)*108)]
% inds{2}=[floor(rand(5,1)*100)+15]
% inds{3}=[floor(rand(5,1)*100)+50]
% col{1}='k'
% col{2}=[.4 .4 1]
% col{3}=[.5 .5 .5]
% for ii=1:length(muinds)
%     plot(pitchtms,pitchdata{ii}(:,inds{ii}),'Color',col{ii})
%     hold on
% end
% % plot(pitchtms,pitchdata{3}(:,58),'c')
% 
% %plot two vertical lines at .04 and .056
% x1=.04
% y1=.056
% plot([.04 .04], [2150 2650],'k--')
% plot([.056 .056], [2150 2650],'k--')
% 
% 
% 
% 
% %generate postdata