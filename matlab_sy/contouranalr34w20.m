%this script is written to generate pitch contours for bk61w42.
clear fv fbins bt

% avls.pvls{1}='/oriole6/bk57w35/test9/'
% avls.cvl{1}='batch.notcatch'
% avls.NT{1}='a'

% tvls.pvls{1}='/oriole6/bk57w35/617_125ma_100Hz_40ms/'
tvls.pvls{1}='/oriole6/r34w20/stimtest724_14bilat_70ma_200Hz_40ms/'
tvls.cvl{1}='batch'
tvls.NT{1}='a'

tvls.pvls{2}='/oriole6/r34w20/14left_100ma_40ms_200Hz/'
tvls.cvl{2}='batch'
tvls.NT{2}='a'

tvls.pvls{3}='/oriole6/r34w20/14_right_80ma_200Hz_40ms_90msdel/'
tvls.cvl{3}='batch'
tvls.NT{3}='a'

tvls.pvls{4}='/oriole6/r34w20/stim724_left23_right14_80ma_40ms_200Hz/'
tvls.cvl{4}='batch'
tvls.NT{4}='a'

tvls.pvls{5}='/oriole6/r34w20/bilat_12_200Hz_80ma_40ms'
tvls.cvl{5}='batch'
tvls.NT{5}='b'

tvls.pvls{6}='/oriole6/r34w20/bilat_34_200Hz_80ma_40ms'
tvls.cvl{6}='batch'
tvls.NT{6}='b'

tvls.pvls{7}='/oriole6/r34w20/bilat_14_run2_100ma_40Hz_200Hz_90msdel'
tvls.cvl{7}='batch'
tvls.NT{7}='b'


tvls.pvls{8}='/oriole6/r34w20/CORRECT2_stim725_bilat14_90msdel_60ms_200Hz_60ma'
tvls.cvl{8}='batch'
tvls.NT{8}='b'

tvls.pvls{9}='/oriole6/r34w20/stim725_100ma_bilat14_200Hz_40ms'
tvls.cvl{9}='batch'
tvls.NT{9}='b'

tvls.pvls{9}='/oriole6/r34w20/test725_newevtaf'
tvls.cvl{9}='batchlim'
tvls.NT{9}='b'

tvls.pvls{10}='/oriole6/r34w20/wnon726'
tvls.cvl{10}='batch.catch'
tvls.NT{10}='b'


tvls.pvls{11}='/oriole6/r34w20/test725_newevtaf'
tvls.cvl{11}='batchlim'
tvls.NT{11}='b'


tvls.pvls{12}='/oriole6/r34w20/wnon726'
tvls.cvl{12}='batch.catchpre'
tvls.NT{12}='b'


tvls.pvls{13}='/oriole6/r34w20/wnon726'
tvls.cvl{13}='batch.catch'
tvls.NT{13}='b'

tvls.pvls{13}='/oriole6/r34w20/stim729'
tvls.cvl{13}='batch'
tvls.NT{13}='b'

tvls.pvls{14}='/oriole6/r34w20/stim730_90ma'
tvls.cvl{14}='batch'
tvls.NT{14}='b'

muinds=[13]


for ii=1:length(muinds)
    crind=muinds(ii)
    pathvl=tvls.pvls{crind}
    cmd=['cd ' pathvl]
    eval(cmd);
%    load sumdata.mat
    

    tbinshft=-.02;
    NFFT=8192;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    
    
      
            fbins{crind}=[2100 3000];
      
    
    save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
    load BINS_B
    NT=tvls.NT{crind};PRENT='';PSTNT='';
    bt{crind}=tvls.cvl{crind}
    
    fv=findwnote4(bt{crind},NT,PRENT,PSTNT,tbinshft,fbins{crind},8192,1,'obs0');
    pitchdata{crind}=jc_pitchcontourFV(fv,1024,1020,1,fbins{crind}(1),fbins{crind}(2),[ 2 ],'obs0')
    contours=pitchdata{crind}
    save pitchdata.mat contours
 end
initsamptime=512/32000;
initsamptimediff=1/8000;
pitchtms=initsamptime:initsamptimediff:(8192-512)/32000
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