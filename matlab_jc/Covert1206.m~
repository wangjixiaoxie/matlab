load /cardinal5/Covert040810.mat

% Compare FFapv to LearningCTL to see if acute APV offsets are correlated with lasting changes
% Compare LearningCTL to zero to make sure all of song isn't shifting in the learning direction 
postday2=postday;
postday2(15,:)=[3 24];
postday2(19,:)=[1 93];
postday2(21,:)=[1 20];
postday2(22,:)=[1 13];
postday2(23,:)=[1 103];
ind2=[1:2 7 10 12:24];
for j=1:length(ind2)
    i=ind2(j);
    coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
    window=[Experiment(i).onCTL:Experiment(i).offCTL];
    FFbefore(i)=mean(mean(Experiment(i).pitchACpreCTL(window,:)));
    FFafter(i)=mean(mean(Experiment(i).pitchACpostCTL(window,postday2(i,1):postday2(i,2))));
    %% for ones where postday is not the same as for
    %% the Targeted note
    Acute(i)=(mean(mean(Experiment(i).pitchAPVCTL(window,round(end/2):end)))-FFbefore(i));
    Permanent(i)=(FFafter(i)-FFbefore(i));    
    LearningCTL(i)=coef*(FFafter(i)-FFbefore(i));
end

% No relationship between Acute APV-dep FF offset in CTL note and Permanent FF offset in CTL note

    corrcoef(Acute(ind2),Permanent(ind2)) % R = 0.0852
    figure;plot(Acute(ind2),Permanent(ind2),'.','Markersize',15)
    hold on;plot([-80 80],[0 0])
    hold on;plot([0 0],[-80 80])
% CTL notes don't shift in the learning direction (control for general shift of song FF)
corrcoef(LearningCTL(ind2),Learning(ind2)) % R = 0.0038
    figure;plot(LearningCTL(ind2),'.','Markersize',15)
    hold on;plot([-80 80],[0 0])
    hold on;plot([0 0],[-80 80])

%%%%%%%
%% APV variability (CV)
%%%%%%%%
        for j=1:length(ind)
                 i=ind(j);
                 window3=[Experiment(i).on:Experiment(i).off];
                % This is just because some parts are noisy - doesn't
                % affect outcome.  Var is reduced everywhere
                 if i==16;window3=[250:350];end
                 if i==4;window3=[180:230];end
                 if i==7;window3=[250:350];end
                 if i==10;window3=[300:400];end
                 a=mean(std(Experiment(i).pitchACpre(window3,:)'));
                 b=mean(mean(Experiment(i).pitchACpre(window3,:)));
                 VarPre(i)=a/b;
                 a=mean(std(Experiment(i).pitchAPV(window3,end-20:end)'));
                 b=mean(mean(Experiment(i).pitchAPV(window3,end-20:end)));                 
                 VarAPV(i)=a/b;
                 a=mean(std(Experiment(i).pitchACpost(window3,postday(i,1):postday(i,1)+100)'));
                 b=mean(mean(Experiment(i).pitchACpost(window3,postday(i,1):postday(i,1)+100)));
                 VarPost(i)=a/b;
                 if i==19
                    a=mean(std(Experiment(i).pitchACpreCTL(300:400,:)'));
                    b=mean(mean(Experiment(i).pitchACpreCTL(300:400,:)));
                    VarPre(i)=a/b;
                    a=mean(std(pitchAPV(300:400,round(end-20):end)'))
                    b=mean(mean(pitchAPV(300:400,round(end-20):end)));
                    VarAPV(i)=a/b;
                    a=mean(std(Experiment(i).pitchACpostCTL(300:400,postday(i,1):postday(i,1)+100)'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(300:400,postday(i,1):postday(i,1)+100)));
                    VarPost(i)=a/b;
                 end
             end
%%%%%%%%%
% Summary data
%%%%%%%%%
            base=[];
            ind=[1 2 4 5 6 7 8 10 12:24];
            for i=ind
                window=round(median(Experiment(i).TargetingWN));
                EtimeAPV=Experiment(i).timeAPVwn;
                [Esorted,sortedAPV]=sort(EtimeAPV);
                mpitchpre=mean(mean(Experiment(i).pitchACpre(window,:)));
                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                a=Experiment(i).pitchAPV(window,round(end/2):end);
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                c=mpitchpre;
                base(i).data=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(1:round(end/3)))))
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                first(i)=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                b=mean(mean(Experiment(i).pitchACpre(window,:)));
                preAPV(i)=coef*(a-b);
                preAPVnorm(i)=preAPV(i)/c;
                a=mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(end/3):round(2*end/3)))));
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                second(i)=coef*(a-b)/c;
                %third(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(2*end/3):round(3*end/3)))))...
                   % -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                %fourth(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(3*end/4):end))))...
                %    -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                a=mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(end-15:end))));
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                finalpoint(i)=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchACpost(window,1:20)));
                b=mean(mean(Experiment(i).pitchACpre(window,:)));
                post1(i)=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchACpost(window,1:postday(i,2)-1)));
                b=mean(mean(Experiment(i).pitchACpre(window,:)));
                post2(i)=coef*(a-b)/c;
                LearningNorm(ind)=Learning(ind)/c;
            end

% preAPVnorm --- % acute offset in adaptive direction - proportional difference
% finalpoint --- % adaptive change in FF during APV - proportional difference
% post1 --- % adaptive change in FF by the next morning - proportional difference




% EXAMPLE 7
load /cardinal6/CovertAnalysis/Exp7.mat
            window1=[140:170];
            % THR ~= 2200Hz
            figure;hold on;subplot(211);hold on;
            plot(Experiment7.timeACpre(114:end),mean(Experiment7.pitchACpre(window1,114:end)),'k.','Markersize',15)
            ratpre=runningaverage(Experiment7.timeACpre(114:end),40);
            rafpre=runningaverage(mean(Experiment7.pitchACpre(window1,114:end)),40);
            plot(ratpre,rafpre,'Color','g','Linewidth',6)
            plot(Experiment7.timeAPV,mean(Experiment7.pitchAPV(window1,:)),'b.','Markersize',15)
            plot(Experiment7.timeAPVwn,mean(Experiment7.pitchAPVwn(window1,:)),'r.','Markersize',15)
            ratapvwn=runningaverage(Experiment7.timeAPVwn,40);
            rafapvwn=runningaverage(mean(Experiment7.pitchAPVwn(window1,:)),40);
            plot(ratapvwn,rafapvwn,'g','Linewidth',6)
            ratapv=runningaverage(Experiment7.timeAPV,40);
            rafapv=runningaverage(mean(Experiment7.pitchAPV(window1,:)),40);
            plot(ratapv,rafapv,'g','Linewidth',6)
            plot(Experiment7.timeACpost,(mean(Experiment7.pitchACpost(window1,:))),'k.','Markersize',15)
            plot(Experiment7.timeACpost(215:end),(mean(Experiment7.pitchACpost(window1,215:end))),'k.','Markersize',15)
            ratpost=runningaverage(Experiment7.timeACpost(215:end),40);
            rafpost=runningaverage(mean(Experiment7.pitchACpost(window1,215:end)),40);
            ratpost=runningaverage(Experiment7.timeACpost,40);
            rafpost=runningaverage(mean(Experiment7.pitchACpost(window1,:)),40);
            plot(ratpost,rafpost,'g','Linewidth',6)
            %avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment7.pitchACpre(window1,:)));
            plot([4975 5010],[mean(rafpre) mean(rafpre)],'k')
            plot([min(ratapvwn) max(ratapvwn)],[2200 2200],'k')

            window1=[180:340];
            subplot(212);hold on;
            plot(Experiment7.timeACpreCTL(114:end),mean(Experiment7.pitchACpreCTL(window1,114:end)),'k.','Markersize',15)
            ratpre=runningaverage(Experiment7.timeACpreCTL(114:end),40);
            rafpre=runningaverage(mean(Experiment7.pitchACpreCTL(window1,114:end)),40);
            plot(ratpre,rafpre,'Color','g','Linewidth',6)
            plot(Experiment7.timeAPVCTL,mean(Experiment7.pitchAPVCTL(window1,:)),'b.','Markersize',15)
            plot(Experiment7.timeAPVwnCTL,mean(Experiment7.pitchAPVwnCTL(window1,:)),'r.','Markersize',15)
            ratapvwn=runningaverage(Experiment7.timeAPVwnCTL,40);
            rafapvwn=runningaverage(mean(Experiment7.pitchAPVwnCTL(window1,:)),40);
            plot(ratapvwn,rafapvwn,'g','Linewidth',6)
            ratapv=runningaverage(Experiment7.timeAPVCTL,40);
            rafapv=runningaverage(mean(Experiment7.pitchAPVCTL(window1,:)),40);
            plot(ratapv,rafapv,'g','Linewidth',6)
            plot(Experiment7.timeACpostCTL(215:end),(mean(Experiment7.pitchACpostCTL(window1,215:end))),'k.','Markersize',15)
            %plot(Experiment7.timeACpostCTL,(mean(Experiment7.pitchACpostCTL(window1,:))),'k.','Markersize',15)
            ratpost=runningaverage(Experiment7.timeACpostCTL(215:end),40);
            rafpost=runningaverage(mean(Experiment7.pitchACpostCTL(window1,215:end)),40);
            plot(ratpost,rafpost,'g','Linewidth',6)
            %avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment7.pitchACpre(window1,:)));
            plot([4975 5010],[mean(mean(Experiment7.pitchACpreCTL(window1,114:end))) mean(mean(Experiment7.pitchACpreCTL(window1,114:end)))],'k')

%%%%%%%%%% REVISED --- 1.20.11            
            % When is the post point to look at?
            clear VarPre VarAPV VarPost1 VarPost2 DayoneCTL DaytwoCTL
            for i=ind2
                % Get the times
                window1=Experiment(i).onCTL:Experiment(i).offCTL;
%                 if isempty(Experiment(i).timeAPVwnCTL)
%                     EtimeACpost=Experiment(i).timeACpostCTL-min(Experiment(i).timeACpostCTL);  % none of these are a problem
%                 end
            EtimeACpost=Experiment(i).timeACpostCTL-max(Experiment(i).timeAPVwnCTL);  % none of these are a problem
                [Esorted,sortedAC]=sort(EtimeACpost);
                dayindices=find(diff(Esorted)>8);
                % Which ones are day one?
                if Esorted(1)<8
                    if isempty(dayindices)
                        indexedACpost=1:length(sortedAC);
                    else
                        indexedACpost=1:dayindices;
                    end
                    dayone=sortedAC(indexedACpost);
                    if dayone(end)+50<length(sortedAC)
                        daytwo=dayone(end)+1:dayone(end)+50;
                    else
                        daytwo=dayone(end)+1:length(sortedAC);
                    end
                else
                    dayone=[];
                    daytwo=1:50;
                end
                DayoneCTL(i).pitch=mean(Experiment(i).pitchACpostCTL(window1,dayone));
                DayoneCTL(i).time=Esorted(dayone);
                DaytwoCTL(i).pitch=mean(Experiment(i).pitchACpostCTL(window1,daytwo));
                DaytwoCTL(i).time=Esorted(daytwo);
                % variability
%                 if length(mean(Experiment(i).pitchACpreCTL(window1,:)))>50
%                 a=mean(std(Experiment(i).pitchACpreCTL(window1,end-50:end)'));
%                 b=mean(mean(Experiment(i).pitchACpreCTL(window1,end-50:end)));
%                 else
                 a=mean(std(Experiment(i).pitchACpreCTL(window1,:)'));
                b=mean(mean(Experiment(i).pitchACpreCTL(window1,:)));         
%                 end
                VarPre(i)=a/b;
                a=mean(std(Experiment(i).pitchAPVCTL(window1,end-20:end)'));
                b=mean(mean(Experiment(i).pitchAPVCTL(window1,end-20:end)));
                VarAPV(i)=a/b;
                if length(daytwo)>49
                    a=mean(std(Experiment(i).pitchACpostCTL(window1,daytwo(1:50))'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(window1,daytwo(1:50))'));
                else
                    a=mean(std(Experiment(i).pitchACpostCTL(window1,daytwo)'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(window1,daytwo)'));
                end
                VarPost2(i)=a/b;
                if length(dayone)>49
                    a=mean(std(Experiment(i).pitchACpostCTL(window1,dayone(end-49:end))'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(window1,dayone(end-49:end))'));
                    FFafter1=mean(mean(Experiment(i).pitchACpostCTL(window1,dayone(end-49:end))));                            
                    VarPost1(i)=a/b;
                    FFpermanent1(i)=FFafter1-FFbefore(i);
                else
                    VarPost1(i)=0;
                end
            end
            
% What is the variability?
mean(VarAPV(ind2)./VarPre(ind2))
bbb=VarPost2(ind2)./VarPre(ind2);
mean(bbb(~isnan(bbb))) % 0.9645
    % Reduction to 73.27% - 26.7% reduction
% What is the offset?            
    ind3=ind2([1:2 4:8 11:12 14:17]); % These are the ones that recovered in time
% mean(Learning(ind3))=31.50+/- 2.023
% 


            clear VarPre VarAPV VarPost1 VarPost2 Dayone Daytwo
            for i=ind
                % Get the times
                window1=Experiment(i).on:Experiment(i).off;
                EtimeACpost=Experiment(i).timeACpost-max(Experiment(i).timeAPVwn);  % time after APVwn off
                [Esorted,sortedAC]=sort(EtimeACpost);
                dayindices=find(diff(Esorted)>8); % find indices for last of the day
                % Which ones are day one?
                if Esorted(1)<8
                    if isempty(dayindices)
                        indexedACpost=1:length(sortedAC);
                    else
                        indexedACpost=1:dayindices(1);
                    end
                    dayone=sortedAC(indexedACpost);
                    if dayone(end)+50<length(sortedAC)
                        daytwo=dayone(end)+1:dayone(end)+50;
                    else
                        daytwo=dayone(end)+1:length(sortedAC);
                    end
                else
                    dayone=[];
                    daytwo=1:50;
                end
                Dayone(i).pitch=mean(Experiment(i).pitchACpost(window1,dayone));
                Dayone(i).time=Esorted(dayone);
                Daytwo(i).pitch=mean(Experiment(i).pitchACpost(window1,daytwo));
                Daytwo(i).time=Esorted(daytwo);
                if length(dayone)>49
                FFafter1targ(i)=mean(mean(Experiment(i).pitchACpost(window1,dayone(end-49:end)))); 
                else 
                 FFafter1targ(i)=mean(mean(Experiment(i).pitchACpost(window1,dayone)));                    
                end
                if length(daytwo)>49
                    FFafter2targ(i)=mean(mean(Experiment(i).pitchACpost(window1,daytwo(1:50)))); 
                else
                    FFafter2targ(i)=mean(mean(Experiment(i).pitchACpost(window1,daytwo))); 
                end
                FFbeforetarg(i)=mean(mean(Experiment(i).pitchACpre(window1,:)));              
            end
  for i=ind
	coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up   
    Learn1(i)=coef*(FFafter1targ(i)-FFbeforetarg(i));
    Learn2(i)=coef*(FFafter2targ(i)-FFbeforetarg(i));
  end
mean(Learn1(ind3(~isnan(Learn1(ind3)))))
mean(Learn2(ind3))