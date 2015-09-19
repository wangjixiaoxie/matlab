%general analysis that goes in and pulls out values for analysis of
%learning

function [sqout]=smsq_anal(ss)

%first find the indices which have stability analysis.

bhvinds=[];
for ii=1:length(ss)
    if(ss(ii).bhvr)
       bhvinds=[bhvinds ii]; 
    end
end

%now for all these indices and not others,
%extract outnotectvls over the stability_window.
indnum=1;

for crind=bhvinds
        
        cmd=['cd ' ss(crind).pth ';load sumdata.mat']
        eval(cmd);
        
        
        ALLNT=avls.ALLSEQNT
        for expnum=avls.LEARNANALIND
            SQNT=avls.SEQTRGNT{expnum}
            sqout(indnum).ssind=crind;
            
            otnt=outnotect;
            bsind=dvl.basind{expnum};
            wnind=dvl.wnind{expnum};
            if(isfield(dvl,'recind'))
                recind=dvl.rec_sl(expnum);
            else
                recind=[];
            end
            if(ss(crind).stab&&(ss(crind).BASFLAG==0))
                 sqout(indnum).BAS=ss(crind).basprob;
            else
                sqout(indnum).BAS=mean(sum(otnt(bsind,SQNT),2)./sum(otnt(bsind,ALLNT),2))
            end
             sqout(indnum).WN=sum(otnt(wnind,SQNT),2)./sum(otnt(wnind,ALLNT),2)
             if(~isempty(recind))
                sqout(indnum).REC=sum(otnt(recind,SQNT),2)./sum(otnt(recind,ALLNT),2)
             else
                 sqout(indnum).REC=NaN;
             end
        indnum=indnum+1;
        end
end
   
    