load /cardinal/Predict30.mat
load /cardinal/FiguresA/CSs73.mat

CS73=CS73final;

%%% 1. How well does my original model do:
            % Code for "Figure3C" from 'Figures3.m'
                        figure;hold on
                        for i=[1:28]
                            [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                            abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                            left=length(abb);
                            abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                            right=length(abb)-left;
                            t=-1*left:1:right-1;
                            hold all
                            plot(t/8,abb,'Linewidth',2)
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
                        hold on;plot(t/8,mnabb(158:1259),'k','Linewidth',6)
                        ylim([0 1.05])
                        xlim([-25 25])
            % Figure 3G - compare predicted vs. actual -
                            % CS predictions - variability with targeting - centered at targ position
                            for i=1:28
                                aax=CS73(i).data;

                                    a(i)=max(aax(Predict(i).onset:Predict(i).offset));

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
                            %%% FINAL PLOT - generates "A1H"
                            j1=max(mnabb(158:1259));
                            j2=max(mnabbT(158:1000));
                            %j3=max(mnabbNT(158:1259));
                            figure;hold on;
                            plot([t/8;t/8],[mnabb(158:1259)/j1+seabb(158:1259)/j1;mnabb(158:1259)/j1-seabb(158:1259)/j1],'color','k')
                            plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',3) % targ imprecision included
                            xlim([-40 40]);ylim([0 1.05])
                            plot(t/8,mnabb(158:1259)/j1)
                            xlim([-25 25])
                          % How much of the shape is explained by the predictions?
                            % 25ms on either side
                                % r=0.9694
                                % r^2=0.9397

                corrcoef(mnabb(500:900),mnabbT(500:900)) % 0.9861

%% 2. Make predictions of new model (PNS)
m=1;
n=1;
trials=1000;
alphas=[0:0.01:0.05];
betas=[0.01:0.02:0.11];
g=0;

for m=1:length(alphas)
    m
    m
    m
    for n=1:length(betas)
        n
        alpha=alphas(m);
        beta=betas(n);
        alpha=0.0005;
        beta=0.0006;
        % Generate predictions for each experiment
        for i=1:28
            PNS(i).Pitch=[];
            PNS(i).Pitch(1,:)=ones(1,size(Predict(i).ResidAC,1));
            for k=2:trials
                % Randomly select FF performance by adding var to current mean
                randomdrawFF=ceil(rand*size(Predict(i).ResidAC,2));
                FFperformance=(1+Predict(i).ResidAC(1:1700,randomdrawFF))'.*PNS(i).Pitch(k-1,1:1700);
                % Does it escape?
                targets=Predict(i).Targeting(Predict(i).Targeting-32<Predict(i).offset & Predict(i).Targeting-32>Predict(i).onset);
                prcval=30+40*isequal(Predict(i).direction,'up');
                Threshold=1+prctile(Predict(i).ResidAC(round(median(targets-32)),:),prcval);
                randomdrawTarg=round(targets(ceil(rand*size(targets))));
                Targtime=[randomdrawTarg-32:randomdrawTarg+32];
                if mean(FFperformance(Targtime))>Threshold
                    escape=1;
                else
                    escape=0;
                end
                % update mean Pitch
                PNS(i).Pitch(k,1:1700)=PNS(i).Pitch(k-1,1:1700)+escape*beta*(FFperformance-PNS(i).Pitch(k-1,1:1700))...
                    -alpha*(PNS(i).Pitch(k-1,1:1700)-PNS(i).Pitch(1,1:1700));
                % sensitivity to deviation as a scalar instead of contour
%                 PNS(i).Pitch(k,1:1700)=PNS(i).Pitch(k-1,1:1700)+escape*beta*(FFperformance-PNS(i).Pitch(k-1,1:1700))...
%                     -alpha*ones(1,1700)*mean(PNS(i).Pitch(k-1,1:1700)-PNS(i).Pitch(1,1:1700));

                
            end
        end
        %% How well does the new model (PNS) do?
        for i=1:28
            LastFF=mean(PNS(i).Pitch(trials-100:trials,1:1700));
            LastFF=(LastFF-1);
            LastFF=LastFF/max(LastFF(Predict(i).onset+30:Predict(i).offset-30));
            LastPitch(i,1:1700)=LastFF(1:1700);
        end
        clear a
        for i=1:28
            aax=LastPitch(i,:);
            a(i)=max(aax(Predict(i).onset:Predict(i).offset));
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
            abb(i,700-dister1:700)=abs(LastPitch(i,Predict(i).onset:Predict(i).onset+b)/a(i));
            abb(i,700:700+dister2)=abs(LastPitch(i,Predict(i).onset+b:Predict(i).offset)/a(i));
        end
        mnabbPNS=zeros(1,1400);
        for i=1:1400
            ind=find(abb(:,i)>0);
            if ~isempty(ind)
                mnabbPNS(i)=mean(abb(ind,i));
            end
        end
        %mnabbPNS=mnabbPNS/max(mnabbPNS(500:900));
        error(m,n)=mean(abs(mnabbPNS(500:900)-mnabb(500:900)));
        holder1=corrcoef(mnabb(500:900),mnabbPNS(500:900));
        corr(m,n)=holder1(2);
        g=g+1;
        guess(g,:)=mnabbPNS;
    end
end
