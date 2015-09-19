load /bulbul5/TWmusc/bs.mat
clearvars -except bs
bsinds=[1 1 2 2 3 3 3 3 3 6 7 7 7];
ACvls=[0 7;0 7;0 6;12 19;0 10;0 10;0 10;0 10;0 10;0 3;0 8;0 8;0 8];
MUvls=[2 3;5 6;5 6;16 17;2 3;4 5;5 6;6 7;8 9;2 3;2 4;4 6;6 8];
COEF=[ones(1,2) -1 ones(1,10)];
for ii=1:length(bsinds)
    [sumbs,shsall,sumshs, shsrev] = shiftanal13(bs, bsinds(ii));
    eval(['cd ' bs(bsinds(ii)).path])
    eval(['load ' bs(bsinds(ii)).matfilename])
    sumbs=sumbs(bsinds(ii));
    beginAC=avls.wn(1).on(1)+ACvls(ii,1);
    endAC=avls.wn(1).on(1)+ACvls(ii,2);
    beginMU=avls.wn(1).on(1)+MUvls(ii,1);
    endMU=avls.wn(1).on(1)+MUvls(ii,2);
    allshiftruns=[];
    for k=1:length(sumbs.shiftruns)
        allshiftruns=[allshiftruns;sumbs.shiftruns{k}];
    end
    % ACSF
                tvacsf=[];fvacsf=[];allFFs=[];alltimes=[]; 
                for i=1:length(avls.adjvls{1}) % for each run
                    if isempty(find(i==allshiftruns))    % if ACSF (i.e. not MUSC)
                        if min(avls.adjvls{1}{1,i}(:,1))>beginAC & max(avls.adjvls{1}{1,i}(:,1))<endAC % if it is a wn run (i.e. not baseline or reversion)
                            if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                tvacsf=[tvacsf tplus'];
                                fvacsf=[fvacsf FFplus'];
                                allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                            else
                                lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                daybegins=[1 lightsout+1];
                                dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                for j=1:length(lightsout) % for each day
                                    tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                    FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                end
                            end
                        end
                    end
                end
                mnFF=mean(allFFs);
                % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
                [vl,ind2]=sort(alltimes);
                gaps=find(diff(alltimes(ind2))>0.3);
                alltimesawake=alltimes(ind2);
                for i=1:length(gaps)
                    alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
                end
                p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
                slope=p(1); 
                alltimesort=alltimes(ind2);
                allFFsort=allFFs(ind2);
    % MUSCIMOL
                tvmu=[];
                fvmu=[];
                count=0;
                clear tvalspre tvalspost fvalspre fvalspost
                for i=1:length(avls.adjvls{1}) % for each run
                    if ~isempty(find(i==allshiftruns))    % if MUSC
                        onset=min(avls.adjvls{1}{1,i}(:,1));
                        offset=max(avls.adjvls{1}{1,i}(:,1));
                        if onset>beginMU & offset<endMU % if it is a wn run (i.e. not baseline or reversion)
                            count=count+1;
                            reltimespre=alltimesort-onset;
                            reltimespost=alltimesort-offset;
                            tvalspre{count}=alltimesort(find(reltimespre<0))-onset;  
                            tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                            fvalspre{count}=allFFsort(find(reltimespre<0));
                            fvalspost{count}=allFFsort(find(reltimespost>0));
                            tmon(count)=avls.tmon(i);
                            tmoff(count)=avls.tmoff(i);

                        end
                    end
                end
                vls(ii).tvalspre=tvalspre;
                vls(ii).tvalspost=tvalspost;
                vls(ii).fvalspre=fvalspre;
                vls(ii).fvalspost=fvalspost;
                t1=char(tmon);
                col=find(t1==':');
                vls(ii).tmon=(str2num(t1(1:col(1)-1))+str2num(t1(col(1)+1:col(1)+2))/60)/24;
                t2=char(tmoff);
                col=find(t2==':');
                vls(ii).tmoff=(str2num(t2(1:col(1)-1))+str2num(t2(col(1)+1:col(1)+2))/60)/24;
                vls(ii).timegap=vls(ii).tmoff-vls(ii).tmon;
% post-processing
                for i=1:length(tvalspre)
                    vls(ii).Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                    vls(ii).FFpre(i)=median(fvalspre{i}(find(tvalspre{i}*24>-8)));
                    vls(ii).CVpre(i)=jcstd(fvalspre{i}(find(tvalspre{i}*24>-8)))/mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                    vls(ii).Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                    vls(ii).Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                    vls(ii).FFpost(i)=median(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                    vls(ii).CVpost(i)=jcstd(fvalspost{i}(find(tvalspost{i}*24>-8)))/mean(fvalspost{i}(find(tvalspost{i}*24>-8)));    
                    vls(ii).Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
                end
                vls(ii).slope=slope;
                vls(ii).COEF=COEF(ii);
                vls(ii).mnFF=mnFF;
                
end

for i=1:length(vls)
    ACTUAL(i)=vls(i).COEF*(vls(i).FFpost-vls(i).FFpre)/vls(i).mnFF;
    POSCTL(i)=vls(i).COEF*(vls(i).slope*vls(i).timegap)/vls(i).mnFF;
    CVpre(i)=vls(i).CVpre;
    CVpost(i)=vls(i).CVpost;
end



                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap
                CVpost
                CVpre

