%the point of this file is to fill in psvals for inactivplotcompar.
figure
ax(1)=subplot(131)
ax(2)=subplot(132)
ax(3)=subplot(133)

plot_all=1;
plot_asymp=0;
runct=1;
ps=[];

for ii=1:length(bs)
     cmd=['cd ' bs(ii).path 'datasum']
    eval(cmd);
    cmd=['load ' bs(ii).matfilename]
    eval(cmd);
    
    if(plot_all)
        ntind=bs(ii).ntind;
        for blocind=1:length(avls.shiftind{ntind});
            
            ps(runct).ntind=ntind;
            ps(runct).birdind=ii;
            indtoplot=avls.shiftind{ntind}{blocind};
            if(~isempty(indtoplot))
                ps(runct).indtoplot=indtoplot;
                if(avls.acmaxz(ntind,blocind)<0)
                    ps(runct).runtype='do'
                    
                else
                    ps(runct).runtype='up'
                end
                ps(runct).plotextraline=1;
                ps(runct).plotbaseline=0;
                runct=runct+1;
            end
            %check for controlnote
            if(~isempty(bs(ii).contrind))
                ps(runct).ntind=bs(ii).contrind
                ps(runct).birdind=ii;
                indtoplot=avls.shiftind{ntind}{blocind};
                if(~isempty(indtoplot))
                    ps(runct).indtoplot=indtoplot;
                    ps(runct).runtype='co';
                   ps(runct).plotextraline=1;
                   ps(runct).plotbaseline=0;
                    runct=runct+1;
                end
             end
        end
    end          
end
for ii=1:(runct-1)
ps(ii).ax=ax;
end
inactivplotcompar3(avls,ps,bs);
linkaxes(ax);        
  
                