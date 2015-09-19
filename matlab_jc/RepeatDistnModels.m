load C:\cardinal8\Repeats\Modelresults.mat
%%%%%%%% 
% The code below creates 
% predprobPRE(experiment #, model #).data
    % 1. MMA (markov with adaptation) - Blue
    % 2. MMA + hardminimum - Red
    % 3. Gamma - Green
% actualprobPRE(experiment #).data - BLACK
%%%%%%%%%%
lw=2;
figure;hold on;
for i=1:length(ind1)
    subplot(3,3,i);hold on;
    plot(actualprobPRE(ind1(i)).data,'k','Linewidth',4)
    plot(predprobPRE(ind1(i),1).data,'b','Linewidth',lw)
    plot(predprobPRE(ind1(i),2).data,'r','Linewidth',lw)
    plot(predprobPRE(ind1(i),3).data,'g','Linewidth',lw)    
    xlim([0 15])
end
%%%%%%%%%%%
% predprobLEARNED(experiment #, learning model #).data
% actualLEARNED(experiment #).data - BLACK
    % 1-3: start with MMA
    % 1. Change P - BLUE
    % 2. Change A - RED
    % 3. Change P & A - GREEN
figure;hold on;
for i=1:length(ind1)
    subplot(3,3,i);hold on;
    plot(actualprobPRE(ind1(i)).data,'k','Linewidth',lw)    
    plot(actualprobLEARNED(ind1(i)).data,'k','Linewidth',4)
    plot(predprobLEARNED(ind1(i),1).data,'b','Linewidth',lw)
    plot(predprobLEARNED(ind1(i),2).data,'r','Linewidth',lw)
    plot(predprobLEARNED(ind1(i),3).data,'g','Linewidth',lw)    
    xlim([0 15])
end
% 4-6: start with Gamma
    % 4. Change shape (NUMBER OF SYSTEMS) - BLUE
    % 5. Change scale (AVERAGE RELIABILITY OF EACH SYSTEM) - RED 
    % 6. Change shape & scale - GREEN
figure;hold on;
for i=1:length(ind1)
    subplot(3,3,i);hold on;
    plot(actualprobPRE(ind1(i)).data,'k','Linewidth',lw)    
    plot(actualprobLEARNED(ind1(i)).data,'k','Linewidth',4)
    plot(predprobLEARNED(ind1(i),4).data,'b','Linewidth',lw)
    plot(predprobLEARNED(ind1(i),5).data,'r','Linewidth',lw)
    plot(predprobLEARNED(ind1(i),6).data,'g','Linewidth',lw)    
    xlim([0 15])
end

%%%%%%%%%%
distances(ind1,:);
resampleddistances(ind1).data
Learningdistances(ind1,:);
resampledLearningdistances(ind1).data;
for z=1:length(ind1)
    i=ind1(z);
    ResampledMARGIN(z)=prctile(resampledLearningdistances(i).data,95);
    ResampledMEDIAN(z)=prctile(resampledLearningdistances(i).data,50);    
end
% Adequate predictions of learning should do better than 95th prctile of
% resampling.  
figure;hold on;
subplot(121);hold on;
plot(ResampledMEDIAN,'k')
plot(ResampledMARGIN,'k')
plot((Learningdistances(ind1,1)),'b','Linewidth',lw,'LineStyle','--') % change P
plot((Learningdistances(ind1,2)),'b','Linewidth',lw,'LineStyle',':') % change A
plot((Learningdistances(ind1,3)),'b','Linewidth',lw) % change P & A
subplot(122);hold on;
plot(ResampledMEDIAN,'k')
plot(ResampledMARGIN,'k')
plot((Learningdistances(ind1,4)),'r','Linewidth',lw,'LineStyle','--') % change P
plot((Learningdistances(ind1,5)),'r','Linewidth',lw,'LineStyle',':') % change A
plot((Learningdistances(ind1,6)),'r','Linewidth',lw) % change P & A
%%%%
for z=1:length(ind1)
    i=ind1(z);
    rpdata=Exp(i).rplength{1};
    actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
    clear distBS
    for bs=1:1000
        scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
        scrambleRPL=rpdata(scrambleIND);
        predprob=hist(scrambleRPL,[1:1:20])/sum(hist(scrambleRPL,[1:1:20]));       
        distBS(bs)=max(abs(actualprob-predprob))/max([actualprob predprob]);
    end
    resampleddistances(i).data=distBS;
    %%%%%%%%
    % Learning
    %%%%%%%%
    i=ind1(z);
    if Exp(i).JC
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-24)); % if JC
    else
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-1)); % if Evren        
    end
    actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
    clear distBS
    for bs=1:1000
        scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
        scrambleRPL=rpdata(scrambleIND);
        predprob=hist(scrambleRPL,[1:1:20])/sum(hist(scrambleRPL,[1:1:20]));
        distBS(bs)=max(abs(actualprob-predprob))/max([actualprob predprob]);
    end
    resampledLearningdistances(i).data=distBS; 
end


%%%%%%%%%
figure;hold on;
for z=1:length(ind1)
    z
    i=ind1(z);
    subplot(3,3,z);hold on;
    %%%%%%%%%
    % Step 1: Fit pre data with p + alpha MMA
    %%%%%%%%%
            rpdata=Exp(i).rplength{1};
            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            clear pcontinue predprob
            for p=0:0.01:1
                for a=0:0.01:1
                    for n=1:20
                        predprob(n)=(a^((n-2)*(n-1)*0.5))*(p^(n-1))*(1-(a^(n-1))*p); % MMA
                    end
                    count=count+1;
                    dist(count)=max(abs(actualprob-predprob))/max([actualprob predprob]);
                    pval(count)=p;
                    aval(count)=a;
                end
            end

            [b,ind]=min(dist);
            p=pval(ind);
            a=aval(ind);
            for n=1:20
                predprob(n)=(a^((n-2)*(n-1)*0.5))*(p^(n-1))*(1-(a^(n-1))*p); % MMA
            end

            plot(actualprob,'k','Linewidth',2);

            plot(predprob,'b','Linewidth',2)
            distances(i,1)=min(dist);
actualprobPRE(i).data=actualprob;
predprobPRE(i,1).data=predprob;
    %%%%%%%%%
    % Step 2: add hardminimum
    %%%%%%%%%
            rpdata=Exp(i).rplength{1};
            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            for hardminimum=0:1:10;
                for p=0.01:0.01:1
                    for a=0.01:0.01:1
                        for n=1:21
                            if n<hardminimum+1
                                pcontinue(n)=1.0;
                            else
                                pcontinue(n)=(a^(n-hardminimum-1))*p;
                            end
                        end
                        for n=2:21
                            predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
                        end
                        count=count+1;
                        %dist(count)=sum(abs(actualprob-predprob).^2);
                        dist(count)=max(abs(actualprob-predprob))/max([actualprob predprob]);
                        pval(count)=p;
                        aval(count)=a;
                        hardmin(count)=hardminimum;
                    end
                end
            end
            [b,ind]=min(dist);
            p_pre=pval(ind);
            a_pre=aval(ind);
            hardminimum_pre=hardmin(ind);
            for n=1:21
                if n<hardminimum_pre+1
                    pcontinue(n)=1.0;
                else
                    pcontinue(n)=(a_pre^(n-hardminimum_pre-1))*p_pre;
                end
            end
            for n=2:21
                predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
            end
            plot(predprob,'r','Linewidth',2)
            %         figure;hold on;
            %         plot(actualprob,'b');plot(predprob,'r')
            distances(i,2)=min(dist);

predprobPRE(i,2).data=predprob;           
         %%%%%
         % model 1: P CHANGES: How good can you do modeling learning by keeping Adap and Hardmin fixed and just changing P?
         %%%%%
    if Exp(i).JC
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-24)); % if JC
    else
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-1)); % if Evren        
    end

            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            hardminimum=hardminimum_pre;
                for p=0.01:0.01:1
                    a=a_pre;
                        for n=1:21
                            if n<hardminimum+1
                                pcontinue(n)=1.0;
                            else
                                pcontinue(n)=(a^(n-hardminimum-1))*p;
                            end
                        end
                        for n=2:21
                            predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
                        end
                        count=count+1;
                        %dist(count)=sum(abs(actualprob-predprob).^2);
                        dist(count)=max(abs(actualprob-predprob))/max([actualprob predprob]);
                        pval(count)=p;
                        aval(count)=a;
                        hardmin(count)=hardminimum;
                end
            [b,ind]=min(dist);
            Learningdistances(i,1)=min(dist);
            p=pval(ind);
            a=aval(ind);
            hardminimum=hardmin(ind);
            for n=1:21
                if n<hardminimum+1
                    pcontinue(n)=1.0;
                else
                    pcontinue(n)=(a^(n-hardminimum-1))*p;
                end
            end
            for n=2:21
                predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
            end
actualprobLEARNED(i).data=actualprob;
predprobLEARNED(i,1).data=predprob;
         %%%%%
         %%%%%
         % model 2: A CHANGES: How good can you do modeling learning by keeping prob and Hardmin fixed and just changing Adaptation?
         %%%%%
    if Exp(i).JC
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-24)); % if JC
    else
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-1)); % if Evren        
    end

            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            hardminimum=hardminimum_pre;
                p=p_pre;
                    for a=0.01:0.01:1
                        for n=1:21
                            if n<hardminimum+1
                                pcontinue(n)=1.0;
                            else
                                pcontinue(n)=(a^(n-hardminimum-1))*p;
                            end
                        end
                        for n=2:21
                            predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
                        end
                        count=count+1;
                        %dist(count)=sum(abs(actualprob-predprob).^2);
                        dist(count)=max(abs(actualprob-predprob))/max([actualprob predprob]);
                        pval(count)=p;
                        aval(count)=a;
                        hardmin(count)=hardminimum;
                end
            [b,ind]=min(dist);
            Learningdistances(i,2)=min(dist);
            p=pval(ind);
            a=aval(ind);
            hardminimum=hardmin(ind);
            for n=1:21
                if n<hardminimum+1
                    pcontinue(n)=1.0;
                else
                    pcontinue(n)=(a^(n-hardminimum-1))*p;
                end
            end
            for n=2:21
                predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
            end            
predprobLEARNED(i,2).data=predprob;           
         %%%%%
         %%%%%
         % model 3: A and P CHANGE: Only hardmin is fixed.
         %%%%%
    if Exp(i).JC
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-24)); % if JC
    else
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-1)); % if Evren        
    end
            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            hardminimum=hardminimum_pre;
                for p=0.01:0.01:1
                    for a=0.01:0.01:1
                        for n=1:21
                            if n<hardminimum+1
                                pcontinue(n)=1.0;
                            else
                                pcontinue(n)=(a^(n-hardminimum-1))*p;
                            end
                        end
                        for n=2:21
                            predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
                        end
                        count=count+1;
                        %dist(count)=sum(abs(actualprob-predprob).^2);
                        dist(count)=max(abs(actualprob-predprob))/max([actualprob predprob]);
                        pval(count)=p;
                        aval(count)=a;
                        hardmin(count)=hardminimum;
                    end
                end
            [b,ind]=min(dist);
            Learningdistances(i,3)=min(dist);
              p=pval(ind);
            a=aval(ind);
            hardminimum=hardmin(ind);
            for n=1:21
                if n<hardminimum+1
                    pcontinue(n)=1.0;
                else
                    pcontinue(n)=(a^(n-hardminimum-1))*p;
                end
            end
            for n=2:21
                predprob(n-1)=prod(pcontinue(1:1:n-1))*(1-pcontinue(n)); % prob of arriving - prob of continuing contingent on having arrived
            end
          
         %%%%%
predprobLEARNED(i,3).data=predprob;           
    %%%%%%%%%%%%%%%%
    %%% can gamma distn do better?
    % Step 3: Fit pre data with gamma
            rpdata=Exp(i).rplength{1};
            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            for shapeG=0:0.5:20
                for scaleG=0:0.05:5
                    gammadist=gampdf([1:1:20],shapeG,scaleG);
                    count=count+1;
                    %dist(count)=sum(abs(actualprob-predprob).^2);
                    dist(count)=max(abs(actualprob-gammadist))/max([actualprob gammadist]);
                    kval(count)=shapeG;
                    thetaval(count)=scaleG;
                end
            end
            [b,ind]=min(dist);
            shapeGpre=kval(ind);
            scaleGpre=thetaval(ind);
            gammadist=gampdf([1:1:20],shapeGpre,scaleGpre);
            %         figure;hold on;
            plot(gammadist,'g','Linewidth',2)
            %%%%%
            distances(i,3)=min(dist);
predprobPRE(i,3).data=gammadist;           
         %%%%%
         % model 4: SCALE CHANGES: How good can you do modeling learning by keeping shape fixed and just changing scale?
         %%%%%
    if Exp(i).JC
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-24)); % if JC
    else
        rpdata=Exp(i).rplength{2}(find(Exp(i).time{2}>Exp(i).time{2}(end)-1)); % if Evren        
    end

            actualprob=hist(rpdata,[1:1:20])/sum(hist(rpdata,[1:1:20]));
            clear dist
            count=0;
            shapeG=shapeGpre;
                for scaleG=0:0.05:5
                    gammadist=gampdf([1:1:20],shapeG,scaleG);
                    count=count+1;
                    %dist(count)=sum(abs(actualprob-predprob).^2);
                    dist(count)=max(abs(actualprob-gammadist))/max([actualprob gammadist]);
                    kval(count)=shapeG;
                    thetaval(count)=scaleG;
                end
            [b,ind]=min(dist);
            Learningdistances(i,4)=min(dist);
                shapeG=kval(ind);
            scaleG=thetaval(ind);
            gammadist=gampdf([1:1:20],shapeG,scaleG);
            %         figure;hold on;
            plot(gammadist,'g','Linewidth',2)
        
predprobLEARNED(i,4).data=gammadist;                       
         %%%%%
         % Model 5: SHAPE CHANGES: How good can you do modeling learning by keeping scale fixed and just changing shape?
         %%%%%         
            clear dist
            count=0;
            scaleG=scaleGpre;
            for shapeG=0:0.5:20
                    gammadist=gampdf([1:1:20],shapeG,scaleG);
                    count=count+1;
                    %dist(count)=sum(abs(actualprob-predprob).^2);
                    dist(count)=max(abs(actualprob-gammadist))/max([actualprob gammadist]);
                    kval(count)=shapeG;
                    thetaval(count)=scaleG;
            end
            [b,ind]=min(dist);
          Learningdistances(i,5)=min(dist);
            shapeG=kval(ind);
            scaleG=thetaval(ind);
            gammadist=gampdf([1:1:20],shapeG,scaleG);
            %         figure;hold on;
            plot(gammadist,'g','Linewidth',2)
          
predprobLEARNED(i,5).data=gammadist;                     
          %%%%%
         % Model 6: SHAPE AND SCALE CHANGE: 
         %%%%%         
         clear dist
         count=0;
         for scaleG=0:0.05:5
             for shapeG=0:0.5:20
                 gammadist=gampdf([1:1:20],shapeG,scaleG);
                 count=count+1;
                 %dist(count)=sum(abs(actualprob-predprob).^2);
                 dist(count)=max(abs(actualprob-gammadist))/max([actualprob gammadist]);
                 kval(count)=shapeG;
                 thetaval(count)=scaleG;
             end
         end
         [b,ind]=min(dist);
         Learningdistances(i,6)=min(dist);
              shapeG=kval(ind);
            scaleG=thetaval(ind);
            gammadist=gampdf([1:1:20],shapeG,scaleG);
            %         figure;hold on;
            plot(gammadist,'g','Linewidth',2)
       
 predprobLEARNED(i,6).data=gammadist;                      
end
