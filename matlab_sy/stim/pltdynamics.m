
function [vlout]=pltdynamics(avls)
    cmd=['cd ' avls.baspath]
    eval(cmd);
    load sumdata.mat
    
for ii=1:length(otmt)
    crmt=otmt(ii)
    crdy=unqdys(ii)-unqdys(1);
  
    tmedges=[7.5/24 9.5/24 12.5/24 15.5/24 20/24]
    endtms=[7/24 21/24]
    colorvls={'k' 'r'}
    
    pltvls{1}.vls=crmt.ct_vls
    pltvls{2}.mvl=crmt.stslope
    pltvls{2}.vls=crmt.st_vls
    pltvls{1}.mvl=crmt.ctslope
    pltvls{2}.yint=crmt.stint
    pltvls{1}.yint=crmt.ctind
    
    for jj=1:2
           crmvl=pltvls{jj}.mvl
            cryint=pltvls{jj}.yint
            plot(endtms+crdy, endtms*crmvl+cryint,'Color',colorvls{jj})
            %catch
            if(jj==1)    
            %catchvls
                vlout(ii).catch=[endtms*crmvl+cryint]
            else
            %stimvls
                vlout(ii).stim=[endtms*crmvl+cryint]
            end
            hold on;
        for vlsind=1:length(pltvls{jj}.vls)
            crvls=pltvls{jj}.vls{vlsind}
         
            crtm=crdy+tmedges(vlsind);
            crmd=median(crvls)
            crer=std(crvls)./sqrt(length(crvls));
            plot(crtm,crmd,'o','Color',colorvls{jj})
            hold on;
            plot([crtm crtm],[crmd+crer crmd-crer],'Color',colorvls{jj})
            
        end
    end
end