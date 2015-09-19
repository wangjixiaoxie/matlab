load /cardinal5/Covert040810.mat
Exps2(4).pitchAPVCTL
Exps1(10).onCTL=940;
Exps1(10).offCTL=1100;
Exps2(2).onCTL=250;
Exps2(2).offCTL=350;
Exps2(3).onCTL=250;
Exps2(3).offCTL=350;
for i=1:10
if isempty(Exps1(i).timeAPVwnCTL)
Exps1(i).timeAPVwnCTL=Exps1(i).timeAPVwn;
end
end
for i=1:10
if isempty(Exps2(i).timeAPVwnCTL)
Exps2(i).timeAPVwnCTL=Exps2(i).timeAPVwn;
end
end
for i=1:4
if isempty(Exps3(i).timeAPVwnCTL)
Exps3(i).timeAPVwnCTL=Exps3(i).timeAPVwn;
end
end
Experiment=[Exps1 Exps2 Exps3];


%%%%%%%%%% REVISED --- 1.20.11            
            % When is the post point to look at?
            clear VarPreCTL VarAPVCTL VarPost1CTL VarPost2CTL DayoneCTL DaytwoCTL
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
                    if isempty(dayindices) | length(Esorted)<50
                        daytwo=1:length(sortedAC);
                    else
                        daytwo=sortedAC(1:50);
                    end
                        
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
                mpitchpreNT(i)=b;
%                 end
                VarPreCTL(i)=a/b;
                a=mean(std(Experiment(i).pitchAPVCTL(window1,end-20:end)'));
                b=mean(mean(Experiment(i).pitchAPVCTL(window1,end-20:end)));
                VarAPVCTL(i)=a/b;
                if length(daytwo)>49
                    a=mean(std(Experiment(i).pitchACpostCTL(window1,daytwo(1:50))'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(window1,daytwo(1:50))'));
                else
                    a=mean(std(Experiment(i).pitchACpostCTL(window1,daytwo)'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(window1,daytwo)'));
                end
%                 a=mean(mean(Experiment(i).pitchAPVwnCTL(window1,end-20:end)));
%                 b=mean(mean(Experiment(i).pitchAPVCTL(window1,end-20:end)));
%                 deltaAPVwnCTL(i)=a-b;
                VarPost2CTL(i)=a/b;
                if length(dayone)>49
                    a=mean(std(Experiment(i).pitchACpostCTL(window1,dayone(end-49:end))'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(window1,dayone(end-49:end))'));
                    FFafter1=mean(mean(Experiment(i).pitchACpostCTL(window1,dayone(end-49:end))));                            
                    VarPost1CTL(i)=a/b;
                    FFpermanent1(i)=FFafter1-mpitchpreNT(i);
                else
                    VarPost1CTL(i)=0;
                end
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
               clear VarPre VarAPV VarPost1 VarPost2 Dayone Daytwo
            for i=ind
                % Get the times
                
                window1=Experiment(i).on:Experiment(i).off;
                if i==4;window1=260:310;end % due to outlier
                if i==10;window1=300:400;end
                if i==16;window1=[250:350];end
%                 if isempty(Experiment(i).timeAPVwnCTL)
%                     EtimeACpost=Experiment(i).timeACpostCTL-min(Experiment(i).timeACpostCTL);  % none of these are a problem
%                 end
            EtimeACpost=Experiment(i).timeACpost-max(Experiment(i).timeAPVwn);  % none of these are a problem
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
                Dayone(i).pitch=mean(Experiment(i).pitchACpost(window1,dayone));
                Dayone(i).time=Esorted(dayone);
                Daytwo(i).pitch=mean(Experiment(i).pitchACpost(window1,daytwo));
                Daytwo(i).time=Esorted(daytwo);
                % variability
%                 if length(mean(Experiment(i).pitchACpreCTL(window1,:)))>50
%                 a=mean(std(Experiment(i).pitchACpreCTL(window1,end-50:end)'));
%                 b=mean(mean(Experiment(i).pitchACpreCTL(window1,end-50:end)));
%                 else
                 a=mean(std(Experiment(i).pitchACpre(window1,:)'));
                b=mean(mean(Experiment(i).pitchACpre(window1,:)));         
%                 end
                VarPre(i)=a/b;
                a=mean(std(Experiment(i).pitchAPV(window1,end-20:end)'));
                b=mean(mean(Experiment(i).pitchAPV(window1,end-20:end)));
                VarAPV(i)=a/b;
                if length(daytwo)>49
                    a=mean(std(Experiment(i).pitchACpost(window1,daytwo(1:50))'));
                    b=mean(mean(Experiment(i).pitchACpost(window1,daytwo(1:50))'));
                else
                    a=mean(std(Experiment(i).pitchACpost(window1,daytwo)'));
                    b=mean(mean(Experiment(i).pitchACpost(window1,daytwo)'));
                end
                VarPost2(i)=a/b;
                if length(dayone)>49
                    a=mean(std(Experiment(i).pitchACpost(window1,dayone(end-49:end))'));
                    b=mean(mean(Experiment(i).pitchACpost(window1,dayone(end-49:end))'));
                    FFafter1=mean(mean(Experiment(i).pitchACpost(window1,dayone(end-49:end))));                            
                    VarPost1(i)=a/b;
                else
                    if length(dayone)>20
                    a=mean(std(Experiment(i).pitchACpost(window1,dayone(end-19:end))'));
                    b=mean(mean(Experiment(i).pitchACpost(window1,dayone(end-19:end))'));
                    FFafter1=mean(mean(Experiment(i).pitchACpost(window1,dayone(end-19:end))));                            
                    VarPost1(i)=a/b;
                    else
                        VarPost1(i)=0;
                    end

                end
                if i==19
                    a=mean(std(Experiment(i).pitchACpreCTL(300:400,:)'));
                    b=mean(mean(Experiment(i).pitchACpreCTL(300:400,:)));
                    VarPre(i)=a/b;
                    a=mean(std(pitchAPV(300:400,round(end-20):end)'))
                    b=mean(mean(pitchAPV(300:400,round(end-20):end)));
                    VarAPV(i)=a/b;
                end
            end
         
            baseCTL=[];

            for i=ind2
                window=Experiment(i).onCTL:Experiment(i).offCTL;
                EtimeAPV=Experiment(i).timeAPVCTL;
                [Esorted,sortedAPV]=sort(EtimeAPV);
                mpitchpre=mean(mean(Experiment(i).pitchACpreCTL(window,:)));
                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                a=Experiment(i).pitchAPVCTL(window,round(end/2):end);
                b=mean(mean(Experiment(i).pitchAPVCTL(window,round(end/2):end)));
                c=mpitchpre;
                baseCTL(i).data=coef*(a-b)/c;
            end
            
     % 95th
     CIbase=[];
     for i=ind
            CIbase(i,:)=resampledCI(base(i).data,0.05,1);
     end
     CIbase=mean(CIbase(ind,:));
            
            