% correctforn.m

% Resampling approach

for i=14%indLes
    CVratio(i)=median(std(Predictx(i).ResidINA(Predict(i).onset:Predict(i).offset,:)')./std(Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:)'));
    %nAC(i)=round(CVratio(i)^2*size(Predictx(i).ResidINA,2));
    %nAC(i)=size(Predictx(i).ResidINA,2);
    gamma=0.5;
    nAC(i)=round(size(Predictx(i).ResidINA,2)*(CVratio(i)/(gamma+(1-gamma)*CVratio(i)))^2);
    if isequal(Predict(i).direction,'up')
        prctile=60;
    else
        prctile=40;
    end

    for k=1:20
        sample=ceil(rand(nAC(i),1)*size(Predict(i).ResidAC,2));
        CSguess=mean(ContingSim2(Predict(i).Targeting,Predict(i).ResidAC(:,sample),prctile));
        CSguess=CSguess-mean(Predict(i).ResidAC(:,sample)');
        CSguesses(k).data=CSguess;
        % normalize and calculate distance
        if i<24 % long notes
                            j=i
                                count=0;
                                for kk=0.5:0.01:2
                                    count=count+1;
                                    choices(count)=sum(abs((kk*abs(CSguess(Predict(i).onset:Predict(i).onset+350))./max(abs(CSguess(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                                end
                                dister(i,k)=min(choices);
                            
        else % short notes
                            j=i
                                count=0;
                                for kk=0.5:0.01:2
                                    count=count+1;
                                    choices(count)=sum(abs((kk*abs(CSguess(Predict(i).onset:Predict(i).onset+160))./max(abs(CSguess(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                                end
                                dister(i,k)=min(choices);
        end
    end
end


Normscore(1:22,1:20)=(dister(1:22,1:20)/350);
Normscore(23:28,1:20)=(dister(23:28,1:20)/160);
for i=1:28
    ScoreN(i)=median(Normscore(i,:));
end

for i=indLes
        if isequal(Predict(i).direction,'up')
        prctile=60;
         else
        prctile=40;
        end
        CSguessLMAN=mean(ContingSim2(Predict(i).Targeting,Predictx(i).ResidINA,prctile));
        if i<24 % long notes
                            j=i
                                count=0;
                                for kk=0.5:0.01:2
                                    count=count+1;
                                    choices(count)=sum(abs((kk*abs(CSguessLMAN(Predict(i).onset:Predict(i).onset+350))./max(abs(CSguessLMAN(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                                end
                                disterLMAN(i)=min(choices);
                            
        else % short notes
                            j=i
                                count=0;
                                for kk=0.5:0.01:2
                                    count=count+1;
                                    choices(count)=sum(abs((kk*abs(CSguessLMAN(Predict(i).onset:Predict(i).onset+160))./max(abs(CSguessLMAN(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                                end
                                disterLMAN(i)=min(choices);
        end
end
% for i=25
%     if isequal(Predict(i).direction,'up')
%         prctile=60;
%     else
%         prctile=40;
%     end
% 
%     for k=1:20
%         sample=ceil(rand(43,1)*size(Predictx(i).ResidINA,2));
%         CSguess=mean(ContingSim2(Predict(i).Targeting,Predictx(i).ResidINA(:,sample),prctile));
%         CSguess=CSguess-mean(Predictx(i).ResidINA(:,sample)');
%         CSguesses(k).data=CSguess;
%         % normalize and calculate distance
%         if i<24 % long notes
%                             j=i
%                                 count=0;
%                                 for kk=0.5:0.01:2
%                                     count=count+1;
%                                     choices(count)=sum(abs((kk*abs(CSguess(Predict(i).onset:Predict(i).onset+350))./max(abs(CSguess(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
%                                 end
%                                 dister25(i,k)=min(choices);
%                             
%         else % short notes
%                             j=i
%                                 count=0;
%                                 for kk=0.5:0.01:2
%                                     count=count+1;
%                                     choices(count)=sum(abs((kk*abs(CSguess(Predict(i).onset:Predict(i).onset+160))./max(abs(CSguess(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
%                                 end
%                                 dister25(i,k)=min(choices);
%         end
%     end
% end


ScoreLM(1:22)=(disterLMAN(1:22)/350);
ScoreLM(23:28)=(disterLMAN(23:28)/160);
%ScoreLM(25)=median(dister25(25,:))/160;

[h,p]=ttest(ScoreLM(indLes)-ScoreN(indLes))
figure;plot(ScoreLM(indLes)-ScoreN(indLes),'*')


%%%%
for k=1:9
    i=indLes(k)
    clear ina
    clear ac
    for j=1:size(Predictx(i).ResidINA,2)
        ina(j,:)=(xcorr(Predictx(i).ResidINA(Predict(i).onset:Predict(i).onset+350,j)));
        ina(j,:)=ina(j,:)/max(ina(j,:));
    end
    for j=1:size(Predictx(i).ResidAC,2)
        ac(j,:)=(xcorr(Predictx(i).ResidAC(Predict(i).onset:Predict(i).onset+350,j)));
        ac(j,:)=ac(j,:)/max(ac(j,:));
    end
    plot(mean(ac))
    hold on;plot(mean(ina),'r')
    %distbt(k)=mean(mean(ac)-mean(ina));
    % fwhm
    distbt(k)=(350-min(find(mean(ac)>0.5)))/(350-min(find(mean(ina)>0.5)));
end
for k=10:13
    i=indLes(k)
    clear ina
    clear ac
    for j=1:size(Predictx(i).ResidINA,2)
        ina(j,:)=(xcorr(Predictx(i).ResidINA(Predict(i).onset:Predict(i).onset+160,j)));
        ina(j,:)=ina(j,:)/max(ina(j,:));
    end        
    for j=1:size(Predictx(i).ResidAC,2)
        ac(j,:)=(xcorr(Predictx(i).ResidAC(Predict(i).onset:Predict(i).onset+160,j)));
        ac(j,:)=ac(j,:)/max(ac(j,:));
    end
    plot(mean(ac))
    hold on;plot(mean(ina),'r')
    %distbt(k)=mean((mean(ac)./max(mean(ac)))-(mean(ina)./max(mean(ina))));
    %distbt(k)=mean(mean(ac)-mean(ina));
    distbt(k)=(160-min(find(mean(ac)>0.5)))/(160-min(find(mean(ina)>0.5)));
end
figure;plot(distbt*100-100,ScoreLM(indLes)-ScoreN(indLes),'*')
%r=0.6
% r^2=0.36
% corrcoef p (non-direction) is 0.03

figure;hold on;subplot(211);plot(1-CVratio(indLes),distbt*100-100,'*')
hold on;plot(1-CVratio(indLes(10:13)),distbt(10:13)*100-100,'*','Color','r')
subplot(212);plot(1-CVratio(indLes),ScoreLM(indLes)-ScoreN(indLes),'*')
hold on;plot(1-CVratio(indLes(10:13)),ScoreLM(indLes(10:13))-ScoreN(indLes(10:13)),'*','Color','r')

%%% Crossco method
        % long notes
            for ii=1:43
                    for k=1:9
                        i=indLes(k);
                        first=ii*8;
                        middle=Predictx(i).onset+first;
                        for j=1:350
                            ab=corrcoef(Predictx(i).ResidAC(j+Predictx(i).onset,:),Predictx(i).ResidAC(middle,:));
                            crosscoAC(j,k)=ab(2);
                            ai=corrcoef(Predictx(i).ResidINA(j+Predictx(i).onset,:),Predictx(i).ResidINA(middle,:));
                            crosscoINA(j,k)=ai(2);
                        end
                        init=350-first+1;
                    end
                    ccAC(ii,init:init+349)=mean(crosscoAC(:,1:9)');
                    ccINA(ii,init:init+349)=mean(crosscoINA(:,1:9)');
                    sdAC(ii,init:init+349)=std(crosscoAC(:,1:9)')/3;
                    sdINA(ii,init:init+349)=std(crosscoINA(:,1:9)')/3;
            end
            for i=1:size(ccAC,2)
                ind1=find(ccAC(:,i)>0);
                mnccAC(i)=mean(ccAC(ind1,i));
                sdccAC(i)=mean(sdAC(ind1,i));
                ind2=find(ccINA(:,i)>0);
                mnccINA(i)=mean(ccINA(ind2,i));
                sdccINA(i)=mean(sdINA(ind1,i));
            end
        % Short notes - ignore for now
            clear crosscoAC
            clear crosscoINA
            for ii=1:20
                    for k=10:13
                        i=indLes(k);
                        first=ii*8;
                        middle=Predictx(i).onset+first;
                        for j=1:160
                            ab=corrcoef(Predictx(i).ResidAC(j+Predictx(i).onset,:),Predictx(i).ResidAC(middle,:));
                            crosscoAC(j,k)=ab(2);
                            ai=corrcoef(Predictx(i).ResidINA(j+Predictx(i).onset,:),Predictx(i).ResidINA(middle,:));
                            crosscoINA(j,k)=ai(2);
                        end

                    end
                        init=350-first+1;
                        ccAC2(ii,init:init+159)=mean(crosscoAC(:,10:13)');
                        ccINA2(ii,init:init+159)=mean(crosscoINA(:,10:13)');
            end
            for i=1:size(ccAC2,2)
                ind1=find(ccAC2(:,i)>0);
                mnccAC2(i)=mean(ccAC2(ind1,i));
                ind2=find(ccINA2(:,i)>0);
                mnccINA2(i)=mean(ccINA2(ind2,i));
            end
            mnBothAC=mnccAC.^2;
            stdBothAC=std(ccAC.^2/sqrt(9));
            mnBothAC(200:500)=mean([mnccAC(200:500).^2;mnccAC2(200:500).^2]);
            stdBothAC(200:500)=std([ccAC(200:500).^2;ccAC2(200:500).^2])/sqrt(13);
            mnBothINA=mnccINA.^2;
            stdBothINA=std(ccINA.^2/sqrt(9));
            mnBothINA(200:500)=mean([mnccINA(200:500).^2;mnccINA2(200:500).^2]);
            stdBothINA(200:500)=std([ccINA(200:500).^2;ccINA2(200:500).^2])/sqrt(13);
%             figure;plot(mnBothAC.^2,'k','Linewidth',2);hold on;plot(mnBothINA.^2,'r','Linewidth',2)
%             xlim([200 499])
%             ylim([0 1])
%             figure;hold on;
%             subplot(121);hold on;
%             plot(mnccAC.^2,'k','Linewidth',2);plot(mnccINA.^2,'r','Linewidth',2)
%             xlim([0 700]);ylim([0 1])
%             subplot(122);hold on;
%             plot(mnccAC2.^2,'k','Linewidth',2);plot(mnccINA2.^2,'r','Linewidth',2)
%             xlim([0 700]);ylim([0 1])

  % Figure 9D
        t1=[1:1:size(mnBothAC,2)];
        figure;hold on;
        plot([t1;t1],[(mnBothAC-stdBothAC);(mnBothAC+stdBothAC)],'-','Color','r')
        plot([t1;t1],[(mnBothINA-stdBothINA);(mnBothINA+stdBothINA)],'-','Color','b')
        xlim([200 499]);ylim([0 1])
  % Figure 7B
        t1=[1:1:size(mnccAC,2)];
        figure;hold on;
        plot([t1;t1],[(mnccINA-sdccINA).^2;(mnccINA+sdccINA).^2],'-','Color','b')
        plot([t1;t1],[(mnccAC-sdccAC).^2;(mnccAC+sdccAC).^2],'-','Color','r')
        xlim([7 692]);ylim([0 1])

% individual variation
clear crosscoAC
clear crosscoINA
clear ccAC
clear ccINA
clear mnccAC
clear mnccINA
            for k=1:9
                    for ii=1:43
                        i=indLes(k);
                        first=ii*8;
                        middle=Predictx(i).onset+40+first;
                        init=350-first+1;
                        for j=1:350
                            ab=corrcoef(Predictx(i).ResidAC(j+Predictx(i).onset+40,:),Predictx(i).ResidAC(middle,:));
                            crosscoAC(init+j,ii)=ab(2);
                            ai=corrcoef(Predictx(i).ResidINA(j+Predictx(i).onset+40,:),Predictx(i).ResidINA(middle,:));
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
            clear crosscoAC
            clear crosscoINA
            for k=10:13
                    for ii=1:18
                        i=indLes(k);
                        first=ii*8;
                        middle=Predictx(i).onset+40+first;
                        init=350-first+1;
                        for j=1:160
                            ab=corrcoef(Predictx(i).ResidAC(j+Predictx(i).onset+40,:),Predictx(i).ResidAC(middle,:));
                            crosscoAC(init+j,ii)=ab(2);
                            ai=corrcoef(Predictx(i).ResidINA(j+Predictx(i).onset+40,:),Predictx(i).ResidINA(middle,:));
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

            for i=1:size(mnccAC,1)
                ind1=find(mnccAC(i,:)>0);
                mnccAC2(i)=mean(mnccAC(i,ind1).^2);
                stdAC2(i)=std(mnccAC(i,ind1).^2);
                ind2=find(mnccINA(i,:)>0);
                mnccINA2(i)=mean(mnccINA(i,ind1).^2);
                stdINA2(i)=std(mnccINA(i,ind1).^2);
            end
            sdAC=stdAC2/sqrt(9);
            sdAC(200:500)=stdAC2(200:500)/sqrt(13);
            sdINA=stdINA2/sqrt(9);
            sdINA(200:500)=stdINA2(200:500)/sqrt(13);
        t1=[1:1:length(mnccAC2)];
        figure;hold on;
        plot([t1;t1],[(mnccAC2-sdAC);(mnccAC2+sdAC)],'-','Color','r')
        plot([t1;t1],[(mnccINA2-sdINA);(mnccINA2+sdINA)],'-','Color','b')
        xlim([17 685]);ylim([0 1])

for i=1:13
    [b,a]=min(find(mnccAC(:,i).^2>0.5));
    [d,c]=min(find(mnccINA(:,i).^2>0.5));
    widthAC(i)=(350-b)/8;
    widthINA(i)=(350-d)/8;
end
figure;plot([1;2],[widthAC;widthINA],'k')
hold on;plot(1,widthAC,'.','MarkerSize',15,'Color','k')
plot(2,widthINA,'.','MarkerSize',15,'Color','r')
xlim([0.5 2.5])
ylim([0 1])


clear crosscoAC
clear crosscoINA
clear ccAC
clear ccINA
clear mnccAC
clear mnccINA
            for k=1:13
                i=indLes(k);
                notelength=Predictx(i).offset-Predictx(i).onset-40;
                numms=floor(notelength/8);
                    for ii=1:numms
                        first=ii*8;
                        middle=Predictx(i).onset+20+first;
                        init=500-first+1;
                        for j=1:notelength
                            ab=corrcoef(Predictx(i).ResidAC(j+Predictx(i).onset+20,:),Predictx(i).ResidAC(middle,:));
                            crosscoAC(init+j,ii)=ab(2);
                            ai=corrcoef(Predictx(i).ResidINA(j+Predictx(i).onset+20,:),Predictx(i).ResidINA(middle,:));
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



%         for k=10:13
%             i=indLes(k);
%             middle=Predictx(i).onset;
%             for j=1:160
%                 ab=corrcoef(Predictx(i).ResidAC(j+Predictx(i).onset,:),Predictx(i).ResidAC(middle,:));
%                 crosscoAC(j+95,k)=ab(2);
%                 ai=corrcoef(Predictx(i).ResidINA(j+Predictx(i).onset,:),Predictx(i).ResidINA(middle,:));
%                 crosscoINA(j+95,k)=ai(2);
%             end
%         end
% 
%         figure;subplot(211);plot(mean(crosscoAC(:,1:9)'))
%         hold on;plot(mean(crosscoINA(:,1:9)'),'r')
%         plot(mean(crosscoAC(:,1:9)')+std(crosscoAC(:,1:9)')/3)
%         plot(mean(crosscoAC(:,1:9)')-std(crosscoAC(:,1:9)')/3)
%         plot(mean(crosscoINA(:,1:9)')+std(crosscoINA(:,1:9)')/3,'r')
%         plot(mean(crosscoINA(:,1:9)')-std(crosscoINA(:,1:9)')/3,'r')
%         subplot(212);hold on;
%         plot(mean(crosscoAC(:,10:13)'))
%         plot(mean(crosscoINA(:,10:13)'),'r')
%         plot(mean(crosscoAC(:,10:13)')+std(crosscoAC(:,10:13)')/2)
%         plot(mean(crosscoAC(:,10:13)')-std(crosscoAC(:,10:13)')/2)
%         plot(mean(crosscoINA(:,10:13)')+std(crosscoINA(:,10:13)')/2,'r')
%         plot(mean(crosscoINA(:,10:13)')-std(crosscoINA(:,10:13)')/2,'r')
%%% PSD
        for i=1:9
            ons=Predictx(indLes(i)).onset;
            [x,psdvalues]=jcpsd2(Predictx(indLes(i)).ResidAC(ons:ons+350,:),8000);
            psdv(i,:)=median(psdvalues);
            [x,psdvaluesINA]=jcpsd2(Predictx(indLes(i)).ResidINA(ons:ons+350,:),8000);
            psdi(i,:)=median(psdvaluesINA);
            ratioPSD(i,:)=median(psdvalues)./median(psdvaluesINA);
        end
%         figure;plot(x,mean(ratioPSD))
%         hold on;plot(x,mean(ratioPSD)+std(ratioPSD)/3)
%         plot(x,mean(ratioPSD)-std(ratioPSD)/3)
        for i=10:13
            ons=Predictx(indLes(i)).onset;
            [x2,psdvalues]=jcpsd2(Predictx(indLes(i)).ResidAC(ons:ons+160,:),8000);
            psdv2(i,:)=median(psdvalues);
            [x2,psdvaluesINA]=jcpsd2(Predictx(indLes(i)).ResidINA(ons:ons+160,:),8000);
            psdi2(i,:)=median(psdvaluesINA);

            ratioPSD2(i-9,:)=median(psdvalues)./median(psdvaluesINA);
        end
%         plot(x,mean(ratioPSD2),'r')
%         hold on;plot(x,mean(ratioPSD2)+std(ratioPSD2)/2,'r')
%         plot(x,mean(ratioPSD2)-std(ratioPSD2)/3,'r')
        %%%
        %%%
        figure;subplot(221);plot(x,psdv','*','Color','b')
hold on;plot(x,psdi','*','Color','r')
hold on;plot(x,psdi','r')
hold on;plot(x,psdv','b')
plot(x2,psdv2','*','Color','b')
hold on;plot(x2,psdi2','*','Color','r')
hold on;plot(x2,psdi2','r')
hold on;plot(x2,psdv2','b')

%%%%%%
        subplot(222);
        hold on;plot(x2,ratioPSD2','Color','r')
        plot(x,ratioPSD','Color','b')
        plot(x2,ratioPSD2','*','Color','r')
        plot(x,ratioPSD','*','Color','b')
      %% GROUPING / STAT ANALYSIS  
            figure;plot(1,[ratioPSD(:,1);ratioPSD(:,2);ratioPSD2(:,1)],'*')
            hold on;plot(2,[ratioPSD(:,3);ratioPSD2(:,2)],'*')
            [h,p]=ttest2([ratioPSD(:,3);ratioPSD2(:,2)],[ratioPSD(:,5);ratioPSD2(:,3)])
            [h,p]=ttest2([ratioPSD(:,3);ratioPSD2(:,2)],[ratioPSD(:,5);ratioPSD(:,4);ratioPSD2(:,3)])
            hold on;plot(3,[ratioPSD(:,5);ratioPSD2(:,3);ratioPSD(:,4)],'*')
            [h,p]=ttest2([ratioPSD(:,1);ratioPSD(:,2);ratioPSD2(:,1)],[ratioPSD(:,3);ratioPSD2(:,2)])
        