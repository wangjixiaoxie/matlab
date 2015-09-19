%wrbnitten 9.20.08
%called by mumeta_compar, and by mure_compar

%principal input is a struct
% fields are ps(run_num).birdind
%            ps(run_num).indtoplot
%            ps(run_num).runtype
%            ps(run_num).ntind
             %ps.ax
             %ps.plotextraline

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
       xvls=ct(ctvl)+xspc:xspc:ct(ctvl)+length(acz)/xspc;
       yvl=muz-acz;
       
       [psout]=getstyle2('do');
       
       for jj=1:length(acz)
           if (ps.plotextraline)
                bot1=[yvl(jj)-avls.mustderrz(ntind,indplt) yvl(jj)+avls.mustderrz(ntind,indplt)]
                bot2=[-avls.acerracz(ntind,indplt) avls.acerracz(ntind,indplt)]
                [psextra]=getstyle2('ex');
                psextra(2).pts=[-1 yvl(jj) 1 yvl(jj)];
                psextra(1).pts=[-1 0 1 0]
                psextra(2).fillpts=[-1 bot1(1); 1 bot1(1); 1 bot1(2); -1 bot1(2);]
                psextra(1).fillpts=[-1 bot2(1); 1 bot2(1); 1 bot2(2); -1 bot2(2)]
                ms=[];
                fillcolor(psextra);
                hold on;
                lineplot(psextra,ms);
           end
           
           psout(1).pts=[0.5 0 0.5 yvl(jj)];
           
           lineplot(psout(1),ms);
           hold on;
       end
          
      
       ct(ctvl)=ct(ctvl)+length(acz)/xspc+xgap;
end
       