%this script is written to generate pitch contours for bk61w42.
clear fv fbins bt

% avls.pvls{1}='/oriole6/bk57w35/test9/'
% avls.cvl{1}='batch.notcatch'
% avls.NT{1}='a'

% tvls.pvls{1}='/oriole6/bk57w35/617_125ma_100Hz_40ms/'
tvls.pvls{1}='/oriole6/bk57w35/623_400Hz_midconfig_15msdel_15ms/'
tvls.cvl{1}='batch'
tvls.NT{1}='a'

tvls.pvls{2}='/oriole6/bk57w35/623_400Hz_mid_zerodel_40ms/'
tvls.cvl{2}='batch'
tvls.NT{2}='a'

tvls.pvls{3}='/oriole6/bk57w35/623_400Hz_midconfig_25msdel_15ms/'
tvls.cvl{3}='batch'
tvls.NT{3}='a'

tvls.pvls{4}='/oriole6/bk57w35/603_400Hz_60ms_origconfig/'
tvls.cvl{4}='batch'
tvls.NT{4}='a'

tvls.pvls{5}='/oriole6/bk57w35/619_20ms_zerodel_400Hz_100ma'
tvls.cvl{5}='batch'
tvls.NT{5}='a'

tvls.pvls{6}='/oriole6/bk57w35/619_middleconfig400Hz_35msdel_30ms_100ma'
tvls.cvl{6}='batch'
tvls.NT{6}='a'

tvls.pvls{7}='/oriole6/bk57w35/619_middlecon_400Hz_20ms_20msdel'
tvls.cvl{7}='batch'
tvls.NT{7}='a'


tvls.pvls{8}='/oriole6/bk57w35/619_middleconfig_400Hz_zerodel_30ms_100ma'
tvls.cvl{8}='batch'
tvls.NT{8}='a'


muinds=[3]


for ii=1:length(muinds)
    crind=muinds(ii)
    pathvl=tvls.pvls{crind}
    cmd=['cd ' pathvl]
    eval(cmd);
%    load sumdata.mat
    

    tbinshft=0;
    NFFT=8192;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    
    
      
            fbins{crind}=[2100 2600];
      
    
    save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
    load BINS_B
    NT=tvls.NT{crind};PRENT='';PSTNT='';
    bt{crind}=tvls.cvl{crind}
    
    fv=findwnote4(bt{crind},NT,PRENT,PSTNT,tbinshft,fbins{crind},8192,1,'obs0');
    pitchdata{crind}=jc_pitchcontourFV(fv,1024,1020,1,fbins{crind}(1),fbins{crind}(2),[3 ],'obs0')
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