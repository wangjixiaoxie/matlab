
% Normalize all curves to maximal adaptive shift in the B direction
        t=[100:-1:-300];
        clear predicted
        clear actual
            for i=1:8
                WNpos=round(median(DShifts(i).toffset)-32+192);
                coef=1-2*isequal(DShifts(i).dirB,'down');
                absshift=coef*shifted(i).data;
                normshifts(i).data=absshift/max(absshift(DShifts(i).onset:DShifts(i).offset));
                actual(i,:)=normshifts(i).data(WNpos+t);
                %%
                abspredshift=coef*CSpred(i).data;
                prednormshifts(i).data=abspredshift/max(abspredshift(DShifts(i).onset:DShifts(i).offset));
                predicted(i,:)=prednormshifts(i).data(WNpos+t);
            end
% Permutation test to calculate error bars
for j=1:401
    j
    for i=1:1000
        draw=ceil(8*(rand(1,8)));
        permmnACT(i)=mean(actual(draw,j));
        draw=ceil(8*(rand(1,8)));
        permmnPRED(i)=mean(predicted(draw,j));
    end
    perm95actual(j)=prctile(permmnACT,95);
    perm5actual(j)=prctile(permmnACT,5);
    perm95pred(j)=prctile(permmnPRED,95);
    perm5pred(j)=prctile(permmnPRED,5);
end
            figure;hold on;
             scaling=mean(actual(:,101));
            scalingP=mean(predicted(:,101));
            plot([t;t],[perm5actual;perm95actual],'k')  
            plot([t;t],[perm5pred;perm95pred],'r')
          
            scaling=mean(actual(:,101));
            scalingP=mean(predicted(:,101));
            plot([t;t],[mean(actual)/scaling-(std(actual)/sqrt(8))/scaling;mean(actual)/scaling+(std(actual)/sqrt(8))/scaling],'k')  
            plot([t;t],[mean(predicted)/scalingP-(std(predicted)/sqrt(8))/scalingP;mean(predicted)/scalingP+(std(predicted)/sqrt(8))/scalingP],'r')
%             ciplot(mean(actual)/scaling-(std(actual)/sqrt(8))/scaling,mean(actual)/scaling+(std(actual)/sqrt(8))/scaling,t,'k',0.8)
%             ciplot(mean(predicted)/scalingP-(std(predicted)/sqrt(8))/scalingP,mean(predicted)/scalingP+(std(predicted)/sqrt(8))/scalingP,t,'r',0.8)
            % Run code for 3C and 3G
            t2=[-400:1:559];
            plot([t2;t2],[mnabb(300:1259)/mnabb(700)-seabb(300:1259)/mnabb(700);mnabb(300:1259)/mnabb(700)+seabb(300:1259)/mnabb(700)],'g')
            plot([0 0],[-2 1.5],'k')
            plot([-192 -192],[-2 1.5],'k')
            plot([-400 300],[0 0],'k')
            xlim([-300 80])
            figure;hold on;
            plot(t,mean(actual)/scaling,'k','Linewidth',3)
            plot(t2,mnabb(300:1259)/mnabb(700),'g','Linewidth',3)
            plot([0 0],[-2 1.5],'k')
            plot([-192 -192],[-2 1.5],'k')
            plot([-400 300],[0 0],'k')
            width=100;
            decayDouble=actual(:,100+width);%./actual(:,101);
            prepoint=[700-width];
            whichind=abb(:,prepoint)~=0;
            decaySingle=abb(whichind,prepoint);%./abb(whichind,700);
            plot(-width,decayDouble,'V','Markersize',10,'Linewidth',1,'Color','k')
            plot(-width,decaySingle,'V','Markersize',10,'Linewidth',1,'Color','g')
            ylim([-2 1.5])
            xlim([-300 80])
            width=80;
            decayDouble=actual(:,100+width);%./actual(:,101);
            prepoint=[700-width];
            whichind=abb(:,prepoint)~=0;
            decaySingle=abb(whichind,prepoint);%./abb(whichind,700);
            plot(-width,decayDouble,'V','Markersize',10,'Linewidth',1,'Color','k')
            plot(-width,decaySingle,'V','Markersize',10,'Linewidth',1,'Color','g')

% for i=1:8
% plot([-300:1:100]/8,normshifts(i).data(round(median(DShifts(i).toffset)-32+192)-300:round(median(DShifts(i).toffset)-32+192)+100),'r','Linewidth',2)
% end
xlim([-20 5])
% Figure 6C - mean +/-SE
                figure;hold on;
                     actual=zeros(8,2000);
                    predicted=zeros(8,2000);

                    for i=1:length(DShifts)
                        if isequal(DShifts(i).dirB,'up')
                            factor=1;
                        else
                            factor=-1;
                        end
                        shifted(i).data=(mean(DShifts(i).pitchALL(:,end-200:end)')-mean(DShifts(i).pitchBaseline'))
                        ons=-192+(DShifts(i).onset-median(DShifts(i).toffset));
                        x=ons:1:ons+(DShifts(i).offset-DShifts(i).onset);
                        maximum=max(abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset)));
                        newx=x+400;
                        actual(i,round(newx))=factor*(1/maximum)*shifted(i).data(DShifts(i).onset:DShifts(i).offset);
                        predicted(i,round(newx))=factor*CSpred(i).data(DShifts(i).onset:DShifts(i).offset);
                        predicted(i,round(newx))=predicted(i,round(newx))/max(abs((predicted(i,round(newx)))));
                    end
                    %hold on;plot([-400:1:1599],bjf,'b')
                    figure;hold on;
                    plot([0-32 0-32],[-1.1 1.1],'k')
                    plot([-192-32 -192-32],[-1.1 1.1],'k')
                    plot([-400 300],[0 0],'k')
                    xlim([-330 120])
                    ylim([-1.1 1.1])
                    for i=1:2000
                        mnactual(i)=0;
                        mnpredicted(i)=0;
                        ind=find(actual(:,i)~=0)
                        mnactual(i)=mean(actual(ind,i));
                        mnSE1actual(i)=mean(actual(ind,i))+std(actual(ind,i))/sqrt(length(ind));
                        mnSE2actual(i)=mean(actual(ind,i))-std(actual(ind,i))/sqrt(length(ind));
                        mnpredicted(i)=mean(predicted(ind,i));
                        mnSE1predicted(i)=mean(predicted(ind,i))+std(predicted(ind,i))/sqrt(length(ind));
                        mnSE2predicted(i)=mean(predicted(ind,i))-std(predicted(ind,i))/sqrt(length(ind));

                    end
                    t1=[-400:1:1599];
                    plot([t1;t1],[mnSE2actual;mnSE1actual],'k')  
                    plot([t1;t1],[mnSE2predicted;mnSE1predicted],'r')
        % Run code for 3C and 3G
                t2=[-400:1:559];
                plot([t2;t2],[mnabb(300:1259)/j1-seabb(300:1259)/j1;mnabb(300:1259)/j1+seabb(300:1259)/j1],'b')

%%%%%%%%%5                
% Figure 3C - All learning curves centered around median targeting position

            figure;hold on
            for i=[1:28]
                [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                left=length(abb);
                abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                right=length(abb)-left;
                t=-1*left:1:right-1;
                plot(t/8,abb,'Linewidth',2,'Color','b')
%                 if i==30 || i==31
%                     plot(t/8,abb,'Linewidth',2,'Color','g')
%                 end
            end
            for i=1:28
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
                ind=find(abb([1:28],i)>0);
                if ~isempty(ind)
                    mnabb(i)=mean(abb(ind,i));
                    seabb(i)=std(abb(ind,i))/sqrt(length(ind));
                end
            end
            t=-542:1:559;
            hold on;plot(t/8,mnabb(158:1259),'r','Linewidth',3)
            ylim([0 1.05])
            xlim([-25 25])
        
% Figure 3G - compare predicted vs. actual -
% RUN CODE for 3C first!!!!!!!
                % CS predictions - variability with targeting - centered at targ position
                for i=1:28
                    aax=CS73(i).data;

                        a(i)=max(aax(Predict(i).onset:Predict(i).offset));

                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                    abb=aax(Predict(i).onset:Predict(i).onset+btop);
                    left=length(abb);
                    abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                    right=length(abb)-left;% Figure 3D - median targeting position vs. maximal FF change

    % This code allows stdev bars but doesn't really add anything.
            %         figure;hold on;
            %         for i=1:28
            %            [a,tops(i)]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
            %            middle(i)=(median(Predict(i).Targeting)-Predict(i).onset)/8;
            %            middleplus(i)=(median(Predict(i).Targeting)+std(Predict(i).Targeting)-Predict(i).onset)/8;
            %            middleminus(i)=(median(Predict(i).Targeting)-std(Predict(i).Targeting)-Predict(i).onset)/8;
            %            tops(i)=tops(i)/8;
            %            plot([middleminus(i) middle(i) middleplus(i)],[tops(i) tops(i) tops(i)],'-','Color','k')
            %         end
        for i=1:28
           [a,tops(i)]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
           middle(i)=median(Predict(i).Targeting)-Predict(i).onset;
        end
            % get in milliseconds 
        middle=middle/8;
        tops=tops/8;
        figure;hold on;
        plot(middle,tops,'.','MarkerSize',15,'Color','k')
        plot(middleplus,tops,'.','MarkerSize',15)
        p=polyfit(middle,tops,1);
        t=0:0.1:100;
        hold on;plot(t,p(2)+p(1)*t,'b')
        xlim([0 70])


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
                    abb(i,700-dister1:700)=abs(CS73(i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                    abb(i,700:700+dister2)=abs(CS73(i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                end
                mnabbT=zeros(1,1400);
                seabbT=zeros(1,1400);
                for i=1:1400
                    ind=find(abb(:,i)>0);
                    if ~isempty(ind)
                        mnabbT(i)=mean(abb(ind,i));
                        seabbT(i)=std(abb(ind,i))/sqrt(length(ind));
                    end
                end
                t=-542:1:559;
                %%%%%
                
    % permutationation
% run 3C immediately prior to get abb
for j=158:1259
    j
    count=0;
    permmnsingle=[];
    % create sample
    sample=abb(:,j);
    sample=sample(sample>0);
    if sum(sample)==0
        perm95single(j)=0;
        perm5single(j)=0;
    else
        for i=1:1000
            draw=ceil(length(sample)*(rand(1,length(sample))));
            permmnsingle(i)=mean(sample(draw));

        end
        perm95single(j)=prctile(permmnsingle,95);
        perm5single(j)=prctile(permmnsingle,5);
    end
end

                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                figure;hold on;
                plot([t/8;t/8],[perm5single(158:1259)/j1;perm95single(158:1259)/j1],'color','k')
                plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',3) % targ imprecision included
                xlim([-40 40]);ylim([0 1.05])
                plot(t/8,mnabb(158:1259)/j1)
                xlim([-25 25])
              % How much of the shape is explained by the predictions?
                % 25ms on either side
                    % r=0.9694
                    % r^2=0.9397               
                
                