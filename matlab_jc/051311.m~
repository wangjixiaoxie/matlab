figure;hold on;
for i=ind
    subplot(5,6,i);hold on;
    plot(jc_residuals(Experiment(i).pitchACpre(goodwin(i).data,end-30:end)),'k')
    plot(jc_residuals(Experiment(i).pitchAPV(goodwin(i).data,end-30:end)),'r')
end
% edit Figures3.m


%%%
clear all
load /cardinal5/Covert040810.mat
load /cardinal5/goodwin.mat
% Covert Learning APV experiments
Experiment=[Exps1 Exps2 Exps3];
clear mnccAC1 mnccINA1
mnccAC=zeros(length(ind),1000);
for k=1:length(ind)
    goodwin(i).data=goodwin(i).data(20:end-20);
    i=ind(k);
    datAC=jc_residuals(Experiment(i).pitchACpre(goodwin(i).data,end-30:end));
    datAPV=jc_residuals(Experiment(i).pitchAPV(goodwin(i).data,end-30:end));
    middle=round(length(goodwin(i).data)/2);
    notelength=length(goodwin(i).data);
    clear crosscoAC
    clear crosscoINA
    %     for j=1:notelength
    %         ab=corrcoef(datAC(j,:),datAC(middle,:));
    %         crosscoAC(j)=ab(2);
    %         ai=corrcoef(datAPV(j,:),datAPV(middle,:));
    %         crosscoINA(j)=ai(2);
    %     end
    for j=1:size(datAC,2)
        crosscoAC(:,j)=xcorr(datAC(:,j));
    end
    for j=1:size(datAPV,2)
        crosscoINA(:,j)=xcorr(datAPV(:,j));
    end
    middle=round(size(crosscoAC,1)/2);
    mnccAC1(k,1000-middle+1:1000-middle+size(crosscoAC,1))=mean(crosscoAC')/max(mean(crosscoAC'));
    mnccINA1(k,1000-middle+1:1000-middle+size(crosscoAC,1))=mean(crosscoINA')/max(mean(crosscoINA'));
end

clearvars -except mnccAC1 mnccINA1
load /cardinal6/LMANtest824.mat
clear mnccAC2 mnccINA2
for k=1:length(indLes)
    i=indLes(k);
    goodwin2(i).data=Predict(i).onset:Predict(i).offset;
    datAC=Predict(i).ResidAC(goodwin2(i).data,:);
    datAPV=Predict(i).ResidINA(goodwin2(i).data,:);
    middle=round(length(goodwin2(i).data)/2);
    notelength=length(goodwin2(i).data);
    clear crosscoAC
    clear crosscoINA
    %     for j=1:notelength
    %         ab=corrcoef(datAC(j,:),datAC(middle,:));
    %         crosscoAC(j)=ab(2);
    %         ai=corrcoef(datAPV(j,:),datAPV(middle,:));
    %         crosscoINA(j)=ai(2);
    %     end
    for j=1:size(datAC,2)
        crosscoAC(:,j)=xcorr(datAC(:,j));
    end
    for j=1:size(datAPV,2)
        crosscoINA(:,j)=xcorr(datAPV(:,j));
    end
    middle=round(size(crosscoAC,1)/2);
    mnccAC2(k,1000-middle+1:1000-middle+size(crosscoAC,1))=mean(crosscoAC')/max(mean(crosscoAC'));
    mnccINA2(k,1000-middle+1:1000-middle+size(crosscoAC,1))=mean(crosscoINA')/max(mean(crosscoINA'));
end

% FWHM
for i=1:size(mnccAC1,1)
	FWHMac1(i)=(1/8)*2*(1000-min(find(mnccAC1(i,:)>0.5)));
	FWHMina1(i)=(1/8)*2*(1000-min(find(mnccINA1(i,:)>0.5)));    
end
for i=1:size(mnccAC2,1)
	FWHMac2(i)=(1/8)*2*(1000-min(find(mnccAC2(i,:)>0.5)));
	FWHMina2(i)=(1/8)*2*(1000-min(find(mnccINA2(i,:)>0.5)));    
end

%
figure;hold on;
plot(median(mnccAC1),'b')
plot(median(mnccAC2),'k')
plot(median(mnccINA1),'r')
plot(median(mnccINA2),'g')

figure;plot(FWHMac1,FWHMina1,'.')
hold on;plot(FWHMac2,FWHMina2,'r.')
hold on;plot([0 15],[0 15],'k')



%%%%%%%%%5
%%%%%%%%%
figure;hold on;

for k=1:length(Experiment)
    [b,indUNQ]=(unique([Experiment(k).timeAPV]));
    pitches=[mean(Experiment(k).pitchAPV(Experiment(k).on:Experiment(k).off,:))];
    clear pt anovavec
    count=0;
    for i=1:length(indUNQ)-1
        these=indUNQ(i):1:indUNQ(i+1)-1;
        for j=1:min([4 length(these)])
            count=count+1;
            pt(i,j)=mean(pitches(these(j))) ;
            anovavec(count,1)=pt(i,j);
            anovavec(count,2)=j;
        end
    end
    clear mnptx
    clear septx
    for i=1:size(pt,2)
        mnptx(i)=median(pt(find(pt(:,i)>0),i));
        septx(i)=std(pt(find(pt(:,i)>0),i))/sqrt(length(pt(find(pt(:,i)>0),i)));
    end
    ffcorrels(k).mnptx=mnptx;
    ffcorrels(k).septx=septx;
    ffcorrels(k).ptx=anovavec;
    [ffcorrels(k).pval,ffcorrels(k).table,ffcorrels(k).stats]=anovan(anovavec(:,1),{anovavec(:,2)},'display','off');
    subplot(5,5,k);hold on;
    plot(mnptx+septx)
    plot(mnptx-septx)
    xlim([0 6])
end
clear anf
for i=1:length(Experiment);anf(i)=ffcorrels(i).table{2,6};end


%
% ANOVA
%
figure;hold on;

for k=1:length(Experiment)
    [b,indUNQ]=(unique([Experiment(k).timeACpre]));
    pitches=mean(Experiment(k).pitchACpre(Experiment(k).on:Experiment(k).off,:));
    clear pt anovavec
    count=0;
    for i=1:length(indUNQ)-1
        these=indUNQ(i):1:indUNQ(i+1)-1;
        for j=1:min([4 length(these)])
            count=count+1;
            pt(i,j)=mean(pitches(these(j))) ;
            anovavec(count,1)=pt(i,j);
            anovavec(count,2)=j;
        end
    end
    clear mnptx
    clear septx
    for i=1:size(pt,2)
        mnptx(i)=median(pt(find(pt(:,i)>0),i));
        septx(i)=std(pt(find(pt(:,i)>0),i))/sqrt(length(pt(find(pt(:,i)>0),i)));
    end
    ffcorrels(k).mnptx=mnptx;
    ffcorrels(k).septx=septx;
    ffcorrels(k).ptx=anovavec;
    [ffcorrels(k).pval,ffcorrels(k).table,ffcorrels(k).stats]=anovan(anovavec(:,1),{anovavec(:,2)},'display','off');
    subplot(5,5,k);hold on;
    plot(mnptx+septx)
    plot(mnptx-septx)
    xlim([0 6])
end
clear anf2
for i=1:length(Experiment);anf2(i)=ffcorrels(i).table{2,6};end

