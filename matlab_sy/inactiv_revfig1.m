%rewritten after botched attempt from before

function []=inactiv_revfig1(sumbs)

figstoplot=1


pstm.runstoplot{1}=[1 2 3]
pstm.runstoplot{2}=[3 4 5]
pstm.bsnum=[8 5]

if(ismember(1,figstoplot))
    axtm=plotstim_tmcourse(sumbs,pstm);
end

function [axtm]=plotstim_tmcourse(sumbs,ps)
    bsln=length(ps.bsnum);
     figure
    for ii=1:bsln
        crbs=ps.bsnum(ii);
       
        axtmcourse(ii)=subplot(1,bsln,ii);
        adjx=1;
        plot_tmcourse2(sumbs(crbs),ps.runstoplot{ii},axtmcourse(ii),adjx);
    end