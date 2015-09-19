%this script is written to generate pitch contours for bk61w42.
muinds=[16 17 18]

for ii=1:length(muinds)
    crind=muinds(ii);
    pathvl=avls.pvls{crind}
    cmd=['cd ' pathvl]
    eval(cmd);
    bt=avls.cvl{crind};
    

    tbinshft=0;
    NFFT=8192;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    fbins=[2000 4000];
    save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
    load BINS_B
    NT='a';PRENT='';PSTNT='';

    
       % edges=[6000:75:8000];
    
    fv{ii}=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
    vals{ii}=getvals(fv{ii},1,'trig')
    pitchdata{ii}=jc_pitchcontourFV(fv{ii},1024,1020,1,2200,2800,[ 3],'obs0')
end
initsamptime=512/32000;
initsamptimediff=1/8000;
pitchtms=initsamptime:initsamptimediff:(8192-512)/32000
figure
subplot(211)
%get the right chunk of data to make example plot

[sm,sp,t,f]=evsmooth(fv{3}(58).datt,32000,10,512,0.8,2,100,10000);
imagesc(t+.016,f,log(abs(sp)));syn;ylim([0,1e4]);
hold on;
plot(pitchtms,pitchdata{3}(:,58),'c')



subplot(212)

%generate 3 random pre
inds{1}=[10 51 21 33]
inds{2}=[63 57 67 70 54 63 75 66]
inds{3}=[85 60 54 65]
col{1}='k'
col{2}=[.4 .4 1]
col{3}=[.5 .5 .5]
for ii=1:length(muinds)
    plot(pitchtms,pitchdata{ii}(:,inds{ii}),'Color',col{ii})
    hold on
end
plot(pitchtms,pitchdata{3}(:,58),'c')

%plot two vertical lines at .04 and .056
x1=.04
y1=.056
plot([.04 .04], [2150 2650],'k--')
plot([.056 .056], [2150 2650],'k--')




%generate postdata