function [CSs2B,distersSELF,distersMUTUAL]=beast(CSs2,Predict)             
            for j=1:28
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    lengthSELF=Predict(j).offset-Predict(j).onset+1;
                    medtarg=median(NormTargs); % expressed as distance from onset
                    proptarg=medtarg/lengthSELF;
                    distPREself=medtarg;
                    distPOSTself=lengthSELF-medtarg;

                    if isequal('up',Predict(j).direction);value=70;else value=30;end
                    for i=1:28
                        i
                        lengthMUTUAL=Predict(i).offset-Predict(i).onset+1;
                        centerMUTUAL=round(lengthMUTUAL/2)+Predict(i).onset;
                        mutualmedtarg=round(proptarg*lengthMUTUAL);
                        NewTargs=NormTargs-medtarg+mutualmedtarg;
                        distPREmutual=mutualmedtarg;
                        distPOSTmutual=lengthMUTUAL-mutualmedtarg;
                        ind1=find(NewTargs>1);
                        ind2=find(NewTargs(ind1)<1800);
                        %CSs2B(i,j).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidAC,value)));
                        % CALCULATE QUALITY OF PREDICTION
                        ptsbefore=round(min([distPREmutual distPREself]));
                        ptsafter=round(min([distPOSTmutual distPOSTself]));
                        ii=i;
                        jj=j;
%                                     count=0;
%                                     for k=0.5:0.01:2
%                                         count=count+1;

%                                         choicesSELF(count)=sum(abs((k*abs(CSs2(jj,jj).data(round(Predict(jj).onset+medtarg-ptsbefore):round(Predict(jj).onset+medtarg+ptsafter)))./max(abs(CSs2(jj,jj).data(round(Predict(jj).onset+medtarg-ptsbefore):round(Predict(jj).onset+medtarg+ptsafter)))))-Predict(jj).LearnedNorm(round(Predict(jj).onset+medtarg-ptsbefore):round(Predict(jj).onset+medtarg+ptsafter))));
%                                         choicesMUTUAL(count)=sum(abs((k*abs(CSs2B(ii,jj).data(round(Predict(ii).onset+mutualmedtarg-ptsbefore):round(Predict(ii).onset+mutualmedtarg+ptsafter)))./max(abs(CSs2B(ii,jj).data(round(Predict(ii).onset+mutualmedtarg-ptsbefore):round(Predict(ii).onset+mutualmedtarg+ptsafter)))))-Predict(jj).LearnedNorm(round(Predict(jj).onset+medtarg-ptsbefore):round(Predict(jj).onset+medtarg+ptsafter))));
%                                     end
%                                     distersSELF(ii,jj)=min(choicesSELF);
%                                     distersMUTUAL(ii,jj)=min(choicesMUTUAL);
                        numpts(i,j)=ptsafter+ptsbefore+1;
                    end
                end