function [probstay,numstay]=repeatsTOptrans2(Exp) 

        count=0;
            for i=1:length(Exp) % for each bird
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
                        indX=(unistart(end):uniend(end));
                    else if j==2 % second half of last day with wn
                            indX=(unistart(end):uniend(end));
                        else if j==3 % day n+2 or more, where day n is last day with wn
                                indX=(unistart(end):uniend(end));
                            end
                        end
                    end
                    for m=2:20
                        probstay(i,j).rn(m-1)=length(find(Exp(i).rplength{j}(indX)>m-1))/length(find(Exp(i).rplength{j}(indX)>m-2));
                        numstay(i,j).rn(m-1)=length(find(Exp(i).rplength{j}(indX)>m-1));

                    end
                end
            end