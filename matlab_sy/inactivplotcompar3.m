%wrbnitten 9.20.08
%called by mumeta_compar, and by mure_compar

%principal input is a struct
% fields are ps(run_num).birdind
%            ps(run_num).indtoplot
%            ps(run_num).runtype
%            ps(run_num).ntind
             %ps.ax
             %ps.plotextraline
             %ps.plotbaseline

function [ax] = inactivplotcompar2(avls,ps,bs)

xspc=1;
xgap=2;

% ax(1)=subplot(131);
% ax(2)=subplot(132);
% ax(3)=subplot(133);

ct(1:3)=0;
for ii=1:length(ps);
    ax=ps(ii).ax;
    run_nm=ii;
    brdind=ps(run_nm).birdind
    rtype=ps(run_nm).runtype;
    ntind=ps(run_nm).ntind
    indplt=ps(run_nm).indtoplot;
    pltextra=ps(run_nm).plotextraline;
    pltbs=ps(run_nm).plotbaseline
    %get the data
    cmd=['cd ' bs(brdind).path 'datasum']
    eval(cmd);
    cmd=['load ' bs(brdind).matfilename]
    eval(cmd);
   
   %set the axis 
   if (rtype=='do')
       axes(ax(1));
       ctvl=1;
   elseif(rtype=='up')
       axes(ax(3));
       ctvl=3;
   else
       axes(ax(2));
       ctvl=2;
   end

   %plot the lines.  (arrows in illustrator);
   
       acz=avls.acz(ntind,indplt);
       muz=avls.muz(ntind,indplt);
       erracz=avls.acerracz(ntind,indplt);
       errmuz=avls.mustderrz(ntind,indplt);
       xvls=ct(ctvl)+xspc:xspc:ct(ctvl)+length(acz)/xspc;
       yvl=muz-acz;
       
       [psout]=getstyle2('do');
       
       for jj=1:length(acz)
           
           
           psout(1).pts=[xvls(jj) 0 xvls(jj) yvl(jj)];
           ms=[];
           lineplot(psout(1),ms);
           hold on;
            if (pltextra)
               if(yvl(jj)>0) 
                   coeff=1
               else
                   coeff=-1
               end
%                coeff=1;
               %this is the mu shift  
                bot1=[yvl(jj)-coeff*erracz(jj)]
                %this is the acsf shift 
                bot2=[coeff*errmuz(jj)]
                
                [psextra]=getstyle2('ex');
                psextra(2).pts=[xvls(jj)-xspc/4 bot1 xvls(jj)+xspc/4 bot1];
                psextra(1).pts=[xvls(jj)-xspc/4 bot2 xvls(jj)+xspc/4 bot2]
                ms=[];
                lineplot(psextra,ms);
           end
       
       end
        
       if(pltbs)
           plot(xvls,acz,'k--')
       end
       
      
       ct(ctvl)=ct(ctvl)+length(acz)/xspc+xgap;
end
       