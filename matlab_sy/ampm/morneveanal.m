%script starting off with morning/evening analyis
%outputs phsumbs.amz  - list of morning z values for each of the days with
%data in the shift
% %        phsumbs(birdind).am(pvlind).ct_tmrange
%                                     .ctvls
%                                     .tmvls
%                                     


%ps.ct=20;
%ps.amtm=8;
%ps.mutm=19;

function [phsumbs] = morneveanal(phsumbs,bs,bsind,ps)
for ii=1:length(bsind)
    crbsind=bsind(ii);
    [avls,graphvals]=loadsumdata(bs, crbsind);
    crbs=phsumbs(crbsind);
    %modify crbs.allasympruns to inclue acsfruns
    %calcasympallruns
    allasympruns=calcallasympruns(crbs);
    phsumbs(crbsind).allasympruns=allasympruns
%getamvls
for indvl=1:length(crbs.aclist(:,1))
   crindvl=crbs.aclist(indvl,1);
   [ct_tmrange,tmedge,ctvls,tmvls]=getvls(crbs,avls,crindvl,'am',ps);
   phsumbs(crbsind).am(crindvl).ct_tmrange=ct_tmrange;
   phsumbs(crbsind).am(crindvl).tmedge=tmedge;
   phsumbs(crbsind).am(crindvl).ctvls=ctvls;
   phsumbs(crbsind).am(crindvl).tmvl=tmvls;
end

for indvl=1:length(crbs.aclist(:,1));
    if(crbs.aclist(indvl,2)~=0)
      crindvl=crbs.aclist(indvl,2);  
    else
      crindvl=crbs.aclist(indvl,1);
    end
    [ct_tmrange,ctvls,tmvls]=getvls(crbs,avls,crindvl,'pm',ps);
     phsumbs(crbsind).pm(crindvl).ct_tmrange=ct_tmrange;
   phsumbs(crbsind).pm(crindvl).ctvls=ctvls;
   phsumbs(crbsind).pm(crindvl).tmvl=tmvls;
end
   
end 
function [ct_tmrange,tmedge,ctvls,tmvls] = getvls(crbs,avls,crindvl,type,ps)
       crvls=avls.allvls{crbs.ntind}{crindvl}
       if(~isempty(crvls))
           if (type=='am')
                amtm=ps.amtm/24;
                tmind=find(mod(crvls(:,1),1)<amtm);
                tmvls=crvls(tmind,2);
           
                ctvls=crvls(1:ps.ct,2);
                ct_tmrange=crvls(ps.ct,1)-crvls(1,1);
                tmedge=mod((crvls(ps.ct,1)),1)*24;
           else
                pmtm=ps.pmtm/24;
                tmind=find(mod(crvls(:,1),1)>pmtm);
                tmvls=crvls(tmind,2);
                lnvls=length(crvls(:,2));
%                 if(lnvls>=ps.ct)
%                     endvl=ps.ct
%                 else
%                     endvl=lnvls
%                 end
                
                
                if((lnvls-ps.ct)>=1)
                    startvl=lnvls-ps.ct
                else
                    startvl=1;
                end
                ctvls=crvls(startvl:end,2);
                ct_tmrange=crvls(end,1)-crvls(startvl,1)
                tmedge=mod(crvls(startvl,1),1)*24;
           end
       else
           
           ctvls=[];
           ct_tmrange=[];
           tmedge=[];
           tmvls=[];


       end
       
       test=1;
           
    function [allasympruns]=calcallasympruns(crbs)
        asympruns=crbs.asympruns
        if(~isempty(asympruns))
        for ii=1:length(asympruns)
            cr_runs=asympruns{ii}
            if(~isempty(cr_runs))
                allasympruns{ii}=find(crbs.flrtmvec>=crbs.flrtmvec(cr_runs(1))&crbs.flrtmvec<=crbs.flrtmvec(cr_runs(end)));
            else
                allasympruns{ii}=[];
            end
        end
        end
       
           
         
        
    

% 
% function [ctvls,ctrange,tmvls,muvls] = doanal(crbs,ps,type)
%     %first find the am runs
%     %by looking for any aclist(1) runs in the allrun period 
%     for shnm=1:length(crbs.shiftruns)
%         if(type=='shi')
%             %include up to first asymptote run
%             %find first asymptote run.
%             muasyruns=crbs.asyruns{shnm}
%             allruns=phsumbs.allshiftruns{shnm};
%             
%             muindvl=find(muasyruns(1),mulist);
%             lastacvl=avls.aclist(muindvl,2);
%             
%             allrunsind=find(allruns==lastacvl);
%             allrunsanal=allruns(1:allrunsind);
%             
%         elseif (type=='asy')
%             allrunsanal=phsumbs.allasympruns{shnm};
%             
%         end
%         [amruns,pmruns]=getampmruns(allrunsanal)
%         
%    %for each of the am runs
%    [ctvls{1},ctrange{1},
% 
% 
%     end
%         
%         
%         
%         
%     