%pulls out shiftruns from stim sumbs struct.
%plots the values, returns a sumplot struct which has all the plotted
%values as well as bsnum,

%calls plotcompar

function [sumplotout,ax]=plotbarstim(sumbs)


    sumplotout=mksumplotout(sumbs);
    ps=makeps(sumplotout);
    ax=ps.ax;
    plotcompar(ps);
    
 %outputs
 %sumplot.acz
 %sumplot.muz
 %sumplot.bsnum
 %first index is upshift, second index is downshift, third index is
 %baseline.
 function [sumplot]=mksumplotout(sumbs)
     ln=4;
     strvls={'acz' 'muz' 'bsnum' 'acbas' 'mubas'}
     [sumplot]=initsumplot(ln,strvls);
 
 for bsnum=1:length(sumbs)
     bs=sumbs(bsnum);
     
    for runvl=1:length(bs.initshiftind)
        crind=bs.initshiftind(runvl);
        if(bs.acz(crind)>0)
            indnum=1;
        else
            indnum=2;
        end
        sumplot(indnum).acz=[sumplot(indnum).acz bs.acz(crind)];
        sumplot(indnum).muz=[sumplot(indnum).muz bs.muz(crind)];
        sumplot(indnum).bsnum=[sumplot(indnum).bsnum bs.bsnum];
    end
    
    if(isfield(bs,'basruns'))
        if(~isempty(bs.basruns))
            [mnac, mnmu]=calcmean(bs, bs.basruns);
            sumplot(3).acbas=[sumplot(3).acbas mnac];
            sumplot(3).mubas=[sumplot(3).mubas mnmu];
            sumplot(3).bsnum=[sumplot(3).bsnum bs.bsnum];
        end
    end
    end
    
     function [mnac,mnmu]=calcmean(bs,basruns)
         mnac=mean(bs.acz(basruns));
         mnmu=mean(bs.muz(basruns));
         
 


function [sumplot]=initsumplot(ln,strvls)
    for ii=1:ln
        for strind=1:length(strvls)
            crstr=strvls{strind};
            cmd=['sumplot(' num2str(ii) ').' crstr '=[]'];
            eval(cmd);
        end
    end
   
                
function ps=makeps(sumplotout)
    figure
    ps.ax=gca();
    ps.xvecs{1}=ones(length(sumplotout(1).acz),1);
    ps.xvecs{2}=2*ones(length(sumplotout(1).acz),1);
    ps.xvecs{3}=ones*ones(length(sumplotout(2).acz),1);   
    ps.xvecs{4}=2*ones(length(sumplotout(2).acz),1); 
    ps.xvecs{5}=4*ones(length(sumplotout(3).acbas),1);
    ps.xvecs{6}=5*ones(length(sumplotout(3).mubas),1);
    
    ps.yvecs{1}=sumplotout(1).acz;
    ps.yvecs{2}=sumplotout(1).muz;
    ps.yvecs{3}=sumplotout(2).acz;
    ps.yvecs{4}=sumplotout(2).muz;
    ps.yvecs{5}=sumplotout(3).acbas;
    ps.yvecs{6}=sumplotout(3).mubas;
    ps.col='k'
    ps.plotbar=1;
    
