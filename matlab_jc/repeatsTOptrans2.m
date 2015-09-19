function [probstay,numstay,semstay]=repeatsTOptrans2(Exp,i) 

                for j=1:length(Exp(i).averages) % for each set (pre,dur,post)
                    if Exp(i).JC==1
                        unidays=find(diff(Exp(i).time{j})>8); % JC time
                        unistart=[1 unidays+1];
                        uniend=[unidays length(Exp(i).time{j})];
                    else
                        unistart=[];
                        uniend=[];
                        unidays=unique(floor(Exp(i).time{j})); % Ev time
                        for kk=1:length(unidays)
                            unistart(kk)=min(find(floor(Exp(i).time{j})==unidays(kk)));
                            uniend(kk)=max(find(floor(Exp(i).time{j})==unidays(kk)));
                        end
                    end

                    if j==1 % day before wn
                        indX=(1:uniend(end));
                    else if j==2 % last day with wn
                            indX=(unistart(end):uniend(end));
                        else if j==3 % day n+2 or more, where day n is last day with wn
                                indX=(unistart(end):uniend(end));
                            end
                        end
                    end
                    for m=2:20
                        probstay(j).rn(m-1)=length(find(Exp(i).rplength{j}(indX)>m-1))/length(find(Exp(i).rplength{j}(indX)>m-2));
                        numstay(j).rn(m-1)=length(find(Exp(i).rplength{j}(indX)>m-1));
                        %
                        % Create surrogate data sets using resampling
                        clear semstayBS
                        rpdata=Exp(i).rplength{j}(indX);
                        for bs=1:1000
                            scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                            scrambleRPL=rpdata(scrambleIND);
                            ptransBS(bs)=length(find(scrambleRPL>m-1))/length(find(scrambleRPL>m-2));
                        end
                        % The standard deviation of these surrogate data sets is equivalent
                        % to the standard error of the mean for the statistic of
                        % interest
                        semstay(j).rn(m-1)=std(ptransBS);
                    end
                end
