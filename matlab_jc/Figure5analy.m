function Figure5analy(Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff,notesizeshiftms)

%%%% Plot the images --- requires manually going into folders
    [F5bk50.avFrontdown,F5bk50.tFrontdown,F5bk50.fFrontdown]=get_avn('batch.catch','a',0.2,0.5,'','','obs0'); 
    figure;
    subplot(141)
    imagesc(F5bk50.tBase,F5bk50.fBase,log(F5bk50.avBase));syn;ylim([0,1e4]);xlim([0 0.2])
    hold on;
    subplot(142)
    imagesc(F5bk50.tEndup,F5bk50.fEndup,log(F5bk50.avEndup));syn;ylim([0,1e4]);xlim([0 0.2])
    subplot(143)
    imagesc(F5bk50.tFrontup,F5bk50.fFrontup,log(F5bk50.avFrontup));syn;ylim([0,1e4]);xlim([0 0.2])
    subplot(144)
    imagesc(F5bk50.tFrontdown,F5bk50.fFrontdown,log(F5bk50.avFrontdown));syn;ylim([0,1e4]);xlim([0 0.2])
%%%% Plot the pitch curves
    % For bk50
    % Panel 1: avgpitch=mean(Alldata(7).exp(1).selectedpitchcurves');
    % Panel 2: avgpitch=mean(Alldata(7).exp(1).selectedpitchcurves');avgpitch=mean(Alldata(7).exp(17).selectedpitchcurves');
    % Panel 3: avgpitch=mean(Alldata(7).exp(1).selectedpitchcurves');avgpitch=mean(Alldata(9).exp(18).selectedpitchcurves');
    % Panel 4: avgpitch=mean(Alldata(8).exp(1).selectedpitchcurves');avgpitch=mean(Alldata(8).exp(25).selectedpitchcurves');
    for i=1:length(avgpitch)
        xax(i)=((i+128)/8)/1000;
    end
    plot(xax,avgpitch*3,'b')

%%%%%%% Plot the distribution of offsets %%%%%%%%
    [to]=gettoffs(Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff,notesizeshiftms);

    figure;
    toind=find(to(7).dat(1:400)>0);
    [n,xout]=hist(to(7).dat(toind),50);
    for i=1:length(xout)
        xout1(i)=((xout(i)+128)/8)/1000;
    end
    hold on;plot(xout1,(n*200)/max(n)+6800,'g')

    [n,xout]=hist(to(9).dat,50);
    for i=1:length(xout)
        xout1(i)=((xout(i)+128)/8)/1000;
    end
    hold on;plot(xout1,7500-(n*200)/max(n),'g')

    [n,xout]=hist(to(8).dat,50);
    for i=1:length(xout)
        xout1(i)=((xout(i)+128)/8)/1000;
    end
    hold on;plot(xout1,(n*200)/max(n)+6800,'g')

%%%% Zoom in for figure B
    ylim([6000 8000])

