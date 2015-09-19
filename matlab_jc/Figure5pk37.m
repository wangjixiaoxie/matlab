function Figure5pk37

%%%% Plot the images --- requires manually going into folders
    [F5pk37.avBase,F5pk37.tBase,F5pk37.fBase]=get_avn('batch.catch','a',0.2,0.5,'','','obs0'); 
    figure;
    subplot(131)
    imagesc(F5pk37.tBase,F5pk37.fBase,log(F5pk37.avBase));syn;ylim([0,1e4]);xlim([0 0.2])
    hold on;
    [F5pk37.avMidUp,F5pk37.tMidUp,F5pk37.fMidUp]=get_avn('batch.catch','a',0.2,0.5,'','','obs0'); 
    subplot(132)
    imagesc(F5pk37.tMidUp,F5pk37.fMidUp,log(F5pk37.avMidUp));syn;ylim([0,1e4]);xlim([0 0.2])
    subplot(133)
    [F5pk37.avMidDn,F5pk37.tMidDn,F5pk37.fMidDn]=get_avn('batch.catch','a',0.2,0.5,'','','obs0'); 
    imagesc(F5pk37.tMidDn,F5pk37.fMidDn,log(F5pk37.avMidDn));syn;ylim([0,1e4]);xlim([0 0.2])

%%%% Plot the pitch curves
    % For pk37
    avgpitch=mean(Alldata(14).exp(1).selectedpitchcurves');
    avgpitch=mean(Alldata(14).exp(27).selectedpitchcurves'); %avgpitch=mean(Alldata(7).exp(17).selectedpitchcurves');
    avgpitch=mean(Alldata(13).exp(3).selectedpitchcurves'); %avgpitch=mean(Alldata(9).exp(18).selectedpitchcurves');
    hold on;
    for i=1:length(avgpitch)
        xax(i)=((i+128)/8)/1000;
    end
    plot(xax,avgpitch*3,'k')

%%%%%%% Plot the distribution of offsets %%%%%%%%
    [to]=gettoffs(Alldata,1,0.8,2,0);

    figure;hist(to(13).dat)
    
    hold on;
    [n,xout]=hist(to(14).dat(1:150),50);
    for i=1:length(xout)
        xout1(i)=((xout(i)+128)/8)/1000;
    end
    hold on;plot(xout1,(n*400)/max(n)+6800,'g')

    [n,xout]=hist(to(13).dat,50);
    for i=1:length(xout)
        xout1(i)=((xout(i)+128)/8)/1000;
    end
    hold on;plot(xout1,7500-(n*400)/max(n),'g')


%%%% Zoom in for figure B
    ylim([6000 8000])

