%purpose of this fxn is to limit data included in analysis before and after
%muscimol run so that ac value represents a stationary value!
function [tmvc] = timeadj_pretm(avls)
    tmvc=avls.adjtimes;
    aclistpre=avls.aclist(:,1);
    aclistpst=avls.aclist(:,2);
    mulist=avls.mulist;
    
    muruns=find(mulist>0)
    
    for ind=1:length(muruns)
        runvl=muruns(ind);
        muindvl=mulist(runvl);
        acpreind=aclistpre(runvl,1);
        acpstind=avls.aclist(runvl,2);
       
         %find where on the list of triplets this is   
         
              must_tm=tmvc(muindvl,1);
              muen_tm=tmvc(muindvl,2);
              
              
              tmpst_tm=must_tm-avls.pretm;
              if(tmpst_tm>tmvc(acpreind,1))
                tmvc(acpreind,1)=tmpst_tm;
              end
              tmpen_tm=muen_tm+avls.pretm;
              if(tmpen_tm<tmvc(acpstind,2))
                tmvc(acpstind,2)=tmpen_tm;
              end
            
            end
          end
   
   