% Zebra finch waveriness
% load cardinal5/MKZF.mat
% Figure 1 - examples
        % Example 1 - b19o61
        figure;hold on;
        subplot(221);hold on;
        plot([0:1/8:250/8],MKZF(1).pitchPRE(150:400,1:30))
        ylim([2300 2800]); xlim([1 30])
        % subplot(232);hold on;
        % plot([0:1/8:250/8],MKZF(1).pitchDIR(150:400,1:30))
        % ylim([2300 2800]); xlim([1 30])
        subplot(222);hold on;
        plot([0:1/8:250/8],MKZF(1).pitchPOST(150:400,1:30))
        ylim([2300 2800]); xlim([1 30])

        % subplot(323);hold on;
        % plot([0:1/8:300/8],MKZF(2).pitchPRE(650:950,1:30))
        % ylim([1400 1800]); xlim([1 30])
        % subplot(235);hold on;
        % plot([0:1/8:300/8],MKZF(2).pitchDIR(650:950,:))
        % ylim([1400 1800]); xlim([1 30])
        % subplot(324);hold on;
        % plot([0:1/8:300/8],MKZF(2).pitchPOST(650:950,1:30))
        % ylim([1400 1800]); xlim([1 30])
        % 


        % Example 2 - pu92w89
        figure;
        subplot(121);hold on;
        plot([0:1/8:400/8],MKZF(4).pitchPRE(600:1000,1:20))
        ylim([2800 3000]); xlim([1 50])
        % subplot(235);hold on;
        % plot([0:1/8:400/8],MKZF(4).pitchDIR(600:1000,1:30))
        % ylim([2800 3100]); xlim([1 50])
        subplot(122);hold on;
        plot([0:1/8:400/8],MKZF(4).pitchPOST(600:1000,1:20))
        ylim([2800 3000]); xlim([1 50])
% Figure 2 - summary data
% /cardinal/FiguresA/FinalFigures/figure9data.mat
        clear crosscoAC
        clear crosscoINA
        clear ccAC
        clear ccINA
        clear mnccAC
        clear mnccINA
        clear mnccAC2
        clear mnccINA2
        for k=1:11
            i=k;
            notelength=PX(i,2)-PX(i,1);
            numms=floor(notelength/8);
            clear crosscoAC
            clear crosscoINA
            for ii=1:numms
                first=ii*8;
                middle=PX(i,1)+first;
                init=500-first+1;
                for j=1:notelength
                    ab=corrcoef(ZFresidsCTL(i).data(j+PX(i,1),:),ZFresidsCTL(i).data(middle,:));
                    crosscoAC(init+j,ii)=ab(2);
                    ai=corrcoef(ZFresidsINA(i).data(j+PX(i,1),:),ZFresidsINA(i).data(middle,:));
                    crosscoINA(init+j,ii)=ai(2);
                end
            end
            for i=1:size(crosscoAC,1)
                ind1=find(crosscoAC(i,:)>0);
                mnccAC(i,k)=mean(crosscoAC(i,ind1));
                ind2=find(crosscoINA(i,:)>0);
                mnccINA(i,k)=mean(crosscoINA(i,ind2));
            end
        end
        %% THIS TAKES THE MEAN
        for i=1:size(mnccAC,1)
            ind1=find(mnccAC(i,:)>0);
            mnccAC2(i)=mean(mnccAC(i,ind1).^2);
            stdAC2(i)=std(mnccAC(i,ind1).^2);
            ind2=find(mnccINA(i,:)>0);
            mnccINA2(i)=mean(mnccINA(i,ind1).^2);
            stdINA2(i)=std(mnccINA(i,ind1).^2);
        end
        sdAC=stdAC2/sqrt(11);
        sdINA=stdINA2/sqrt(11);
   %% THIS PLOTS IT
        t1=[1:1:length(mnccAC2)];
        figure;hold on;
        plot([t1;t1],[(mnccAC2-sdAC);(mnccAC2+sdAC)],'-','Color','b')
        plot([t1;t1],[(mnccINA2-sdINA);(mnccINA2+sdINA)],'-','Color','r')
        plot([t1;t1],[(mnccINA2-sdINA)*mean(1-CVratio);(mnccINA2+sdINA)*mean(1-CVratio)],'-','Color','k')
        plot([t1;t1],[(mnccAC2-sdAC);(mnccAC2+sdAC)],'-','Color','b')
        xlim([220 780]);ylim([0 1])
        
%%% FIGURES FOR LMAN only component
        for i=1:11
            CVratio(i)=1-mean(std(ZFresidsINA(i).data(PX(i,1):PX(i,2),:)'))/mean(std(ZFresidsCTL(i).data(PX(i,1):PX(i,2),:)'));
        end
        for i=1:size(mnccAC,2)
            weightedFF(:,i)=(mnccAC(:,i).^2-CVratio(i)*mnccINA(:,i).^2)/(1-CVratio(i));
            weightedLesion(:,i)=CVratio(i)*mnccINA(:,i).^2;
            weightedFF(:,i)=weightedFF(:,i)/max(weightedFF(:,i));
        end
        for i=1:size(mnccAC,1)
            ind1=find(weightedFF(i,:)>0);
            weightedFF2(i)=mean(weightedFF(i,ind1));
            MEANLesion2(i)=mean(weightedLesion(i,ind1));
            SELesion2(i)=std(weightedLesion(i,ind1))/sqrt(i);
        end
        weightedFF2=weightedFF2/max(weightedFF2);
        figure;plot(weightedFF2,'Linewidth',2)
        halfwidth=(500-399-min(find(weightedFF2(400:end)>0.5)))/4;
        % 13.215 (precisely)
 %%%%%%%%%5
 %%%%%%%%%%5
 %%%%%%%%%%%
        figure;hold on;
        for i=1:11
            subplot(3,4,i)
            hold on;
            plot([-25:1/8:25],sqrt(mnccINA(300:700,i)),'r','Linewidth',2);xlim([-25 25]);ylim([0 1])
            plot([-25:1/8:25],sqrt(mnccAC(300:700,i)),'b','Linewidth',2);xlim([-25 25]);ylim([0 1])
            plot([-25:1/8:25],sqrt(mnccINA(300:700,i))*(1-CVratio(i)),'k','Linewidth',2);xlim([-25 25]);ylim([0 1])
        end

        figure;hold on;
        mnccAC3=sqrt(mnccAC2);
        mnccINA3=sqrt(mnccINA2);
        sdAC3=(sdAC);
        sdINA3=(sdINA);
        plot([t1;t1],[(mnccAC3-sdAC3);(mnccAC3+sdAC3)],'-','Color','b')
        plot([t1;t1],[(mnccINA3-sdINA3);(mnccINA3+sdINA3)],'-','Color','r')
        plot([t1;t1],[(mnccINA3-sdINA3)*mean(1-CVratio);(mnccINA3+sdINA3)*mean(1-CVratio)],'-','Color','k')
        plot([t1;t1],[(mnccAC3-sdAC3);(mnccAC3+sdAC3)],'-','Color','b')
        xlim([220 780]);ylim([0 1])
%%%%
%%%%
% Sonograms - this looks great!!!
sonoPRE=jc_ifd062310(MKZF(4).rawPRE(1:20,:),1024,1020,2,2600,3100,[1],'obs0',1);
sonoPOST=jc_ifd062310(MKZF(4).rawPOST(1:20,:),1024,1020,2,2600,3100,[1],'obs0',1);
clims=[0 10];
figure;hold on
for i=1
    subplot(2,10,i);hold on;
    imagesc([1/8:1/8:1245/8],[32000:-32000/512:32000/512],log(sonoPRE(i).data),clims);colormap(bone)
    xlim([60 140])
    ylim([2000 8000])
end
for i=2:10
    subplot(2,10,i);hold on;
    set(gca,'xtick',[],'ytick',[])
    imagesc([1/8:1/8:1245/8],[32000:-32000/512:32000/512],log(sonoPRE(i).data),clims);colormap(bone)
    xlim([60 140])
    ylim([2000 8000])
end

for i=1:10
    subplot(2,10,10+i);hold on;
    set(gca,'xtick',[],'ytick',[])
    imagesc([1/8:1/8:1245/8],[32000:-32000/512:32000/512],log(sonoPOST(i).data),clims);colormap(bone)
    xlim([60 140])
    ylim([2000 8000])
end 
%%%%%%%
%%%%%%%5
%%%%%%%%%
% 06.28.10

% op92
        sonoUDpre=jc_ifd062310(UDpre',1024,1020,1,3600,4400,[1],'wav',1);
        sonoUDmusc=jc_ifd062310(UDmusc',1024,1020,1,3600,4400,[1],'wav',1);
        sonoDpre=jc_ifd062310(Dpre',1024,1020,1,3600,4400,[1],'wav',1);
        sonoDmusc=jc_ifd062310(Dmusc',1024,1020,1,3600,4400,[1],'wav',1);
        clims=[0 5];
        figure;hold on
            subplot(1,4,1);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(sonoUDpre(1).data),clims);colormap(bone)
            ylim([0 10000])
            xlim([0 160])
            subplot(1,4,2);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(sonoDpre(1).data),clims);colormap(bone)
            ylim([0 10000])
              xlim([0 160])  
            subplot(1,4,3);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(sonoUDmusc(1).data),clims);colormap(bone)
            ylim([0 10000])
              xlim([0 160])  
            subplot(1,4,4);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(sonoDmusc(1).data),clims);colormap(bone)
            ylim([0 10000])
            xlim([0 160])

        for i=820:1100
            a=corrcoef(pitchpyUDmuAll(i,:),pitchpyUDmuAll(960,:));
            b=corrcoef(pitchpyUDacAll(i,:),pitchpyUDacAll(960,:));
            ccMU(i)=a(2);
            ccAC(i)=b(2);
        end
        for i=200:800
            a=corrcoef(pitchUDmuAll(i,:),pitchUDmuAll(500,:));
            b=corrcoef(pitchUDacAll(i,:),pitchUDacAll(500,:));
            ccMU2(i)=a(2);
            ccAC2(i)=b(2);
        end
        figure;plot(ccMU2,'r')
        hold on;plot(ccAC2,'b')

% py40
        PYsonoUDpre=jc_ifd062310(pyUDacEG',1024,1020,1,3600,4400,[1],'wav',1);
        PYsonoUDmusc=jc_ifd062310(pyUDmuEG',1024,1020,1,3600,4400,[1],'wav',1);
        PYsonoDpre=jc_ifd062310(pyDacEG',1024,1020,1,3600,4400,[1],'wav',1);
        PYsonoDmusc=jc_ifd062310(pyDmuEG',1024,1020,1,3600,4400,[1],'wav',1);
        clims=[0 5];
        figure;hold on
            subplot(1,4,1);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(PYsonoUDpre(1).data),clims);colormap(bone)
            ylim([0 10000])
            xlim([70 150])
            subplot(1,4,2);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(PYsonoDpre(1).data),clims);colormap(bone)
            ylim([0 10000])
              xlim([70 150])  
            subplot(1,4,3);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(PYsonoUDmusc(1).data),clims);colormap(bone)
            ylim([0 10000])
              xlim([70 150])  
            subplot(1,4,4);hold on;
            imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(PYsonoDmusc(1).data),clims);colormap(bone)
            ylim([0 10000])
            xlim([70 150])

        for i=820:1100
            a=corrcoef(pitchpyUDmuAll(i,:),pitchpyUDmuAll(960,:));
            b=corrcoef(pitchpyUDacAll(i,:),pitchpyUDacAll(960,:));
            ccMU(i)=a(2);
            ccAC(i)=b(2);
        end
        for i=200:800
            a=corrcoef(pitchUDmuAll(i,[1:3 5:10]),pitchUDmuAll(500,[1:3 5:10]));
            b=corrcoef(pitchUDacAll(i,[1:3 5:10]),pitchUDacAll(500,[1:3 5:10]));
            ccMU2(i)=a(2);
            ccAC2(i)=b(2);
        end
% Pitch contours        
        figure;hold on;
        subplot(231);hold on;
        for i=1:10
        plot([0:1/8:600/8],pitchUDacAll(200:800,i),'b')
        end
        ylim([4000 4400]); xlim([0 600/8])
        subplot(232);hold on;
        for i=1:10
        plot([0:1/8:600/8],pitchUDmuAll(200:800,i),'r')
        end
        ylim([4000 4400]); xlim([1 600/8])
        subplot(234);hold on;
        for i=1:10
        plot([0:1/8:280/8],pitchpyUDacAll(820:1100,i),'b')
        end
        xlim([0 280/8]); ylim([3250 3650])
        subplot(235);hold on;
        for i=1:10
        plot([0:1/8:280/8],pitchpyUDmuAll(820:1100,i),'r')
        end
        xlim([0 280/8]); ylim([3250 3650])
        subplot(233);hold on;
        plot([0:1/8:600/8],ccMU2(200:800),'k')
        plot([0:1/8:600/8],ccAC2(200:800),'b')
        plot([0:1/8:600/8],ccMU2(200:800)*(mean(std(pitchpyUDmuAll(820:1100,:)))/mean(std(pitchpyUDacAll(820:1100,:)))),'r')
        xlim([0 600/8])
        subplot(236);hold on;
        plot([0:1/8:280/8],ccMU(820:1100),'k')
        plot([0:1/8:280/8],ccAC(820:1100),'b')
        plot([0:1/8:280/8],ccMU(820:1100)*(mean(std(pitchpyUDmuAll(820:1100,:)))/mean(std(pitchpyUDacAll(820:1100,:)))),'r')
        xlim([0 280/8])
%%%%%
figure;hold on
for i=1:10
    subplot(2,10,i);hold on;
    set(gca,'xtick',[],'ytick',[])
    imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(op92sonoUDmusc(i).data),clims);colormap(hot)
    xlim([1 140])
    ylim([1 12000])
end 
for i=1:10
    subplot(2,10,10+i);hold on;
    set(gca,'xtick',[],'ytick',[])
    imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(op92sonoUDacsf(i).data),clims);colormap(hot)
    xlim([1 140])
    ylim([1 12000])
end 






% 7.3.10
% Consensus is 1.3 for all data
op92sonoUDmusc=jc_ifd062310(UDmuAll,1024,1016,1.5,2600,3100,[1],'obs0',1);
op92sonoUDacsf=jc_ifd062310(UDacAll,1024,1016,1.5,2600,3100,[1],'obs0',1);
py40sonoUDmusc=jc_ifd062310(pyUDmuAll,1024,1016,1.5,2600,3100,[1],'obs0',1);
py40sonoUDacsf=jc_ifd062310(pyUDacAll,1024,1016,1.5,2600,3100,[1],'obs0',1);
% py40 - not impressive difference
figure;hold on
for i=1:10
    subplot(2,10,i);hold on;
    set(gca,'xtick',[],'ytick',[])
    imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(py40sonoUDmusc(i).data),clims);colormap(bone)
    xlim([70 150])
    ylim([1 12000])
end 
for i=1:10
    subplot(2,10,10+i);hold on;
    set(gca,'xtick',[],'ytick',[])
    imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(py40sonoUDacsf(i).data),clims);colormap(bone)
    xlim([70 150])
    ylim([1 12000])
end 

%%%% op
clims=[0 5]
figure;hold on;
subplot(2,2,1);
set(gca,'xtick',[],'ytick',[])
imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(op92sonoUDacsf(4).data),clims);colormap(hot)
xlim([20 120])
ylim([1 12500])
syn
subplot(2,2,2);
set(gca,'xtick',[],'ytick',[])
imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(op92sonoUDmusc(1).data),clims);colormap(hot)
syn
xlim([20 120])
ylim([1 12500])
subplot(2,2,3);
set(gca,'xtick',[],'ytick',[])
imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(op92sonoUDacsf(4).data),clims);colormap(hot)
xlim([20 120])
ylim([8000 12500])
syn
subplot(2,2,4);
set(gca,'xtick',[],'ytick',[])
imagesc([1/8:1/8:1245/8],[44100:-44100/512:44100/512],log(op92sonoUDmusc(1).data),clims);colormap(hot)
syn
xlim([20 120])
ylim([8000 12500])


%%% Consensus function
op92medAC=zeros(10,7);
op92medMU=zeros(10,7);
py40medAC=zeros(10,7);
py40medMU=zeros(10,7);
for i=1:10
    i
    op92medAC(i,:)=ConsentPlot2(UDacAll(i,:),44100,1024,1020,1,1,3,3);
    op92medMU(i,:)=ConsentPlot2(UDmuAll(i,:),44100,1024,1020,1,1,3,3);
end
for i=1:10
    i
    py40medAC(i,:)=ConsentPlot2(pyUDacAll(i,:),44100,1024,1020,1,1,3,3);
    py40medMU(i,:)=ConsentPlot2(pyUDmuAll(i,:),44100,1024,1020,1,1,3,3);
end
