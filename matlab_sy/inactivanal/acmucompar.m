function [ratvls]=acmucompar(avls)
    %for every muscimol run, find an acsf pre and post partner (if they
    %exist)
    PRETM=4;
    PSTM=4;
    for ii=1:length(avls.muanal)
        muind=avls.muanal(ii)
        
        strt_tm=avls.rawtimes(muind,1);
        floormu=floor(strt_tm);
        end_tm=avls.rawtimes(muind,2);
        
        
        %find all the ac start timees which are the same day, ordinarily
        %there should be two.
        floorstrt_all=floor(avls.rawtimes(:,1));
        
        
        acpreind=find(floorstrt_all==floormu&avls.rawtimes(:,1)<strt_tm);
        acpstind=find(floorstrt_all==floormu&avls.rawtimes(:,1)>strt_tm);
        
        for ntind=1:length(avls.NT)
            if (~isempty(acpreind))
                %select predata
                ind=find(avls.adjvls{ntind}{acpreind(end)}(:,1)>strt_tm-PRETM/24)
                predata=avls.adjvls{ntind}{acpreind(end)}(ind,:);
            else
                predata=[];
            end
            if(~isempty(acpstind))
                ind=find(avls.adjvls{ntind}{acpstind(1)}(:,1)<end_tm+PSTM/24)
                postdata=avls.adjvls{ntind}{acpstind(1)}(ind,:);
            else
                postdata=[];
            end
                combdata=[predata; postdata]
                           acmn=mean(combdata(:,2))
                           ratvls(muind).acmean{ntind}=acmn
                           acstd=std(combdata(:,2))
                           ratvls(muind).acstd{ntind}=std(combdata(:,2))
%                         ratvls(muind).acstderr{ntind}=std(combdata(:,2))/sqrt(length(combdata(:,2)))
                            mumn=avls.mnvl{ntind}(muind);
                ratvls(muind).mumean{ntind}=avls.mnvl{ntind}(muind);
                mustdv=avls.stdv{ntind}(muind);
                ratvls(muind).mustd{ntind}=avls.stdv{ntind}(muind);
%             ratvls(muind).mustderr{ntind}=avls.stdv{ntind}(muind)/sqrt(length(avls.adjvls{ntind}(muind)));
                cvrat=(mustdv/mumn)/(acstd/acmn);
                ratvls(muind).cvfac=1/cvrat;

        end
    end
        