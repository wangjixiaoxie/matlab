% This agrees pretty well with the no phase method - fakepsd.m

% The one weakness of this method (any method) is that if any component is
% represented more in the LMAN disrupted context (rare), then we can't
% account for that.  

% The method takes the difference between the principle components and
% generates LMAN contour variants based on the CV's (for each component) present in
% the LMAN intact context



% load /cardinal/FiguresA/LMANtest824.mat
%% This loop gets fake data
for k=indLes
    if k>22
        numpts=160;
    else
        numpts=350;
    end
            size1=size(Predictx(k).ResidAC,2);
            size2=size(Predictx(k).ResidINA,2);
            beginning=Predictx(k).onset;
            ending=Predictx(k).onset+numpts;
            width=ending-beginning+1;
            x1=Predictx(k).ResidAC(beginning:ending,:);
            x2=Predictx(k).ResidINA(beginning:ending,:);
            [PC, SCORE, LATENT]=princomp([x1 x2]');
        %     figure;plot(SCORE(:,1),SCORE(:,2),'*')
        %     hold on;plot(SCORE(1:size1,1),SCORE(1:size1,2),'*','Color','r')
                % How much does LMAN contribute to each
                meanAC=mean(abs(SCORE(1:size1,:)));
                meanINA=mean(abs(SCORE(size1+1:end,:)));
                stdAC=std(abs(SCORE(1:size1,:)));
                stdINA=std(abs(SCORE(size1+1:end,:)));
        %     figure;plot((meanAC-meanINA)*PC'/max((meanAC-meanINA)*PC'),'g')
        %     hold on;plot((meanAC*PC')/max(meanAC*PC'),'b')
        %     plot((meanINA*PC')/max(meanINA*PC'),'r')
       
%%% crosscorr method - shows that the fake curves are appropriate
    clear samples
    g=[];
    for j=1:500
        coefs=(abs(randn(1,length(PC)))*(stdAC/meanAC)).*(meanAC-meanINA)+(meanAC-meanINA);
        ind=find(coefs<0);
        coefs(ind)=0;

        for i=1:width
            m=rand;
            if m>0.5
                g(i)=-1;
            else 
                g(i)=1;
            end
        end
        samples(j,:)=coefs.*g*PC';
    end
    LMANcurves=samples';
               notelength=width;
            numms=floor(notelength/8);
            clear crosscoAC
            clear crosscoINA
            clear mnccAC
            clear mnccINA
            clear mnccLMAN
            clear crosscoLMAN
            for ii=1:numms
                first=ii*8;
                firstpoint=beginning;
                middle=firstpoint+first;
                init=500-first+1;
                for j=1:notelength
                    ab=corrcoef(x1(j,:),x1(middle-firstpoint,:));
                    crosscoAC(init+j,ii)=ab(2);
                    ai=corrcoef(x2(j,:),x2(middle-firstpoint,:));
                    crosscoINA(init+j,ii)=ai(2);
                    index1=find(~isnan(mean(samples')));
                    aL=corrcoef(LMANcurves(j,index1),LMANcurves(middle-firstpoint,index1));
                    crosscoLMAN(init+j,ii)=aL(2);
                end
            end
            for i=1:size(crosscoAC,1)
                ind1=find(crosscoAC(i,:)>0);
                mnccAC(i)=mean(crosscoAC(i,ind1));
                ind2=find(crosscoINA(i,:)>0);
                mnccINA(i)=mean(crosscoINA(i,ind2));
                ind3=find(crosscoLMAN(i,:)>0);
                mnccLMAN(i)=mean(crosscoLMAN(i,ind3));
            end

%             figure;plot(mnccAC)
%             hold on;plot(mnccINA,'r')
%             plot(mnccLMAN,'g')
            LMANcomponent(k).data=LMANcurves;
end

for i=indLes
    if i>22
        numpts=160;
    else
        numpts=350;
    end
    
    if isequal(Predictx(i).direction,'up')
        thr=60;
        coef=1;
    else
        thr=40;
        coef=-1;
    end
    Targs=round(Predictx(i).Targeting-(Predictx(i).onset));
    ind1=find(Targs>0);
    ind2=find(Targs(ind1)<numpts);
    CSmock(i).data=ContingSim2(Targs(ind1(ind2)),LMANcomponent(i).data,thr);
    
    %normalize
    count=0;
    predicted=coef*mean(CSmock(i).data);
    for k=0.5:0.01:2
        count=count+1;
        choices(count)=sum(abs((k*abs(predicted)./max(abs(predicted)))-Predictx(i).LearnedNorm(Predictx(i).onset:Predictx(i).onset+numpts)));
    end
    disterLMAN2(i)=min(choices)/numpts;
end
figure;hold on;
plot([1;2],[disterLMAN2(indLes);BestSelfCSnorm64(indLes)],'k')
plot(1,disterLMAN2(indLes),'.','MarkerSize',25,'Color','g')
plot(2,BestSelfCSnorm64(indLes),'.','MarkerSize',25,'Color','r')
ylim([0 0.2]);xlim([0.8 2.2])





%%%%%%%%%
%%%%%%%%
%%%% Actual adaptation for indLes birds
       figure;hold on
            for i=[1:28]
                [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                left=length(abb);
                abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                right=length(abb)-left;
                t=-1*left:1:right-1;
            end
            for i=[1:28]
                notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
            end
            abb=zeros(28,1400);
            for i=[1:28]
                b=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                b=round(b);
                dister1=b;
                dister2=notewidth(i)*8-b;

                abb(i,700-dister1:700)=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+b);
                abb(i,700:700+dister2)=Predict(i).LearnedNorm(Predict(i).onset+b:Predict(i).offset);
            end
            mnabb=zeros(1,1400);
            seabb=zeros(1,1400);
            for i=1:1400
                ind=find(abb(indLes,i)>0);
                if ~isempty(ind)
                    mnabb(i)=mean(abb(indLes(ind),i));
                    seabb(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                end
            end
            t=-542:1:559;
            figure;hold on;plot(t/8,mnabb(158:1259)/max(mnabb(158:1259)),'k','Linewidth',3)
            ylim([0 1.05])
            xlim([-40 40])
%%% Predicted from all data for indLes birds
                for i=1:28
                    aax=CSs64(i,i).data;
                    if isequal(Predict(i).direction,'up')
                        a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                    else
                        a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                    end
                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                    abb=aax(Predict(i).onset:Predict(i).onset+btop);
                    left=length(abb);
                    abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                    right=length(abb)-left;
                    t=-1*left:1:right-1;
                    abb=abb/a(i);
                end
                for i=1:28
                    notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                end
                abb=zeros(28,1400);
                for i=1:28
                    b=round(median(Predict(i).Targeting)-Predict(i).onset); %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                    dister1=(b);
                    dister2=(notewidth(i)*8-b);
                    abb(i,700-dister1:700)=abs(CSs64(i,i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                    abb(i,700:700+dister2)=abs(CSs64(i,i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                end
                mnabbT=zeros(1,1400);
                seabbT=zeros(1,1400);
                for i=1:1400
                    ind=find(abb(indLes,i)>0);
                    if ~isempty(ind)
                        mnabbT(i)=mean(abb(indLes(ind),i));
                        seabbT(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                    end
                end
                t=-542:1:559;
                %%%%%
                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',3) % targ imprecision included
%%% Prediction for just DMP

for i=indLes
    predictionLMAN(i).data=mean(CSmock(i).data);
end
    
%                 for i=indLes
%                     if i>22
%                         numpts=160;
%                     else
%                         numpts=350;
%                     end
%                     aax=predictionLMAN(i).data;
%                     Targs=round(Predictx(i).Targeting-(Predictx(i).onset));
%                     ind1=find(Targs>0);
%                     ind2=find(Targs(ind1)<numpts);
%                     Targs=Targs(ind1(ind2));
% 
%                     if isequal(Predict(i).direction,'up')
%                         a(i)=max(aax);
%                     else
%                         a(i)=min(aax);
%                     end
%                     btop=median(Targs);
%                     abb=aax(1:btop);
%                     left=length(abb);
%                     abb=[abb aax(btop:end)];
%                     right=length(abb)-left;
%                     t=-1*left:1:right-1;
%                     abb=abb/a(i);
%                 end
                for i=indLes
                    notewidth(i)=numpts./8;
                end
                abb=zeros(28,1400);
                for i=indLes
                    if i>22
                        numpts=160;
                    else
                        numpts=350;
                    end

                    Targs=round(Predictx(i).Targeting-(Predictx(i).onset));
                    ind1=find(Targs>0);
                    ind2=find(Targs(ind1)<numpts);
                    Targs=Targs(ind1(ind2));

                    b=round(median(Targs)); %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                    dister1=(b);
                    dister2=(numpts-b);
                    abb(i,700-dister1:700)=abs(predictionLMAN(i).data(1:b+1)/a(i));
                    abb(i,700:700+dister2)=abs(predictionLMAN(i).data(b+1:end)/a(i));
                end
                mnabbT=zeros(1,1400);
                seabbT=zeros(1,1400);
                for i=1:1400
                    ind=find(abb(indLes,i)>0);
                    if ~isempty(ind)
                        mnabbT(i)=mean(abb(indLes(ind),i));
                        seabbT(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                    end
                end
                t=-542:1:559;
                %%%%%
                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                plot(t/8,mnabbT(158:1259)/0.631,'g','LineWidth',3) % targ imprecision included
                
                % DMP only
clear a
                for m=1:13
                    i=indLes(m);
                    aax=CSLMAN(m).data;
                    if isequal(Predict(i).direction,'up')
                        a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                    else
                        a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                    end
                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                    abb=aax(Predict(i).onset:Predict(i).onset+btop);
                    left=length(abb);
                    abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                    right=length(abb)-left;
                    t=-1*left:1:right-1;
                    abb=abb/a(i);
                end
                for i=1:28
                    notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                end
                abb=zeros(28,1400);
                for m=1:13
                    i=indLes(m);
                    b=round(median(Predict(i).Targeting)-Predict(i).onset); %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                    dister1=(b);
                    dister2=(notewidth(i)*8-b);
                    abb(i,700-dister1:700)=abs(CSLMAN(m).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                    abb(i,700:700+dister2)=abs(CSLMAN(m).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                end
                mnabbT=zeros(1,1400);
                seabbT=zeros(1,1400);
                for i=1:1400
                    ind=find(abb(indLes,i)>0);
                    if ~isempty(ind)
                        mnabbT(i)=mean(abb(indLes(ind),i));
                        seabbT(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                    end
                end
                t=-542:1:559;
                %%%%%
                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                plot(t/8,mnabbT(158:1259)/j2,'b','LineWidth',3) % targ imprecision included

