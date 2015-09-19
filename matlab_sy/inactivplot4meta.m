%input parameters are avls,graphvals,muparam, and histogram,axhandle,notevl
%heightoffset
 ax=figure;

    htcnt=0;
    zrat=1;
 for birdind=1:length(bs)      
% length(bs)

    htcnt=htcnt+1.5
    axis([-6000 8000 0 14])
    ntind=bs(birdind).ntind;
    
    cmd=['cd ' bs(birdind).path 'datasum']
    eval(cmd);
    cmd=['load ' bs(birdind).matfilename]
        eval(cmd);
        
    for blocind=1:length(avls.muanal{ntind})
        for dyind=1:length(avls.muanal{ntind}{blocind} )   
            ps.ax=ax
            ps.line=1
            ps.ht(1)=htcnt+1;
            ps.ht(2)=htcnt+1.5;
            ps.hist=[];
            dyindvl=avls.muanal{ntind}{blocind}(dyind)
            ps.zrat=1;
            ps.tickht=.5
            ps.notevl=ntind;
            ps.mu='mu'
            ps.dyind=find(avls.mulist==dyindvl);
            inactivplot4(avls,graphvals,ps);
            htcnt=htcnt+1.5;
        end
    end

text(-5,htcnt,avls.bname)
 end
% plot([0 0],[0 50],'k--');
xlabel('pitch z-score');
box off;