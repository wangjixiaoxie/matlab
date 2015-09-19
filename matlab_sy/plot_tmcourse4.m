%rewritten 11.21.08

% plot_timecourse4.m
%modified 7.29.09
%in order to generalize to stimdata


function [ax] = plot_tmcourse4(sumbs, plotind,STIM)
if ~(exist('STIM'))
    STIM=0
end


for ii=1:length(sumbs)
    ps=setplotvls;
    ps.plotruns=1;
    if(STIM)
        ps.plotruns=0
    end
    ps.plotwn=1;
    ppf=ps.ppf
    plotcount=3;
    
    pn=mod(ii,ppf)
    %if this is the last plot
    if(pn==0)
        pn=ppf;
    end
    if(pn==1)
        figure;
    end
    
    ax(pn)=subplot(ppf,1,pn)
    [vls]=getvalsmod(sumbs(ii),STIM);
    plotvals(vls,ax(pn),ps);
end


    
    function [ps] = setplotvls()
    
        ps.lw=3
        ps.msize=9
        ps.errlw=2
        ps.mucol=[0.4 0.4 1]
        ps.mufillcol=[0.82 0.82 0.82]
        ps.acfillcol=[0.92 0.96 0.98]
        ps.accol='k'
        ps.meanlw=2
        ps.errht=50
        ps.acmark='o'
        ps.mumark='o'
        ps.ppf=3;
  
  function   [vls]=getvalsmod(sumbs);
      vls.ntind=sumbs.ntind;
      ntind=sumbs.ntind;
     
      vls.subvl=floor(sumbs.rawtimes(1,1));
      indvl=find(sumbs.mulist~=0)
      for ii=1:length(indvl)
            ind=indvl(ii)
                if(~isempty(sumbs.mulist(ind)))
                    indvl2=sumbs.mulist(ind);
                    vls.acz(ii)=sumbs.acz(ntind,indvl2);
                    vls.aczerr(ii)=sumbs.acerrz(ntind,indvl2)
                    vls.muz(ii)=sumbs.muz(ntind,indvl2);
                    vls.muzerr(ii)=sumbs.muerrz(ntind,indvl2)
                end
      end
    vls.rwtms=sumbs.rawtimes(sumbs.mulist(indvl),1)
    vls.rawtmsall=sumbs.rawtimes(:,1);
    if(STIM==0)
        vls.basruns=sumbs.basruns;
        vls.revruns=sumbs.revruns;
        vls.shiftruns=sumbs.shiftruns
        vls.initruns=sumbs.initind;
        vls.subruns=sumbs.subruns;
    end
    vls.initmean=sumbs.initmean;
    vls.mulist=sumbs.mulist;
    vls.bname=sumbs.bname;
    vls.wnon=sumbs.flr_wnon
    vls.wnoff=sumbs.flr_wnoff
  function []=plotvals(vls,ax,ps);
    axes(ax);
    rwtm=vls.rwtms
    rwtmall=vls.rawtmsall;
    sbval=vls.subvl;
    acz=vls.acz;
    initmean=vls.initmean;
    bname=vls.bname
    aczer=vls.aczerr
    muz=vls.muz
    muzer=vls.muzerr
    plot(rwtm-sbval,acz,'Color',ps.accol,'Linewidth',ps.lw);
    hold on;
     plot(rwtm-sbval,acz,'Marker',ps.acmark,'Color',ps.accol,'Markersize',ps.msize,'MarkerFaceColor',ps.accol); 
    plot(rwtm-sbval,muz,'Color',ps.mucol,'Linewidth',ps.lw);
    plot(rwtm-sbval,muz,'Marker',ps.mumark,'Color',ps.mucol,'Markersize',ps.msize,'MarkerFaceColor',ps.mucol);
    
    %plot vertical lines as error bars

    plot([(rwtm-sbval)'; (rwtm-sbval)'],[(acz-aczer); (acz+aczer)],'k','Linewidth',ps.errlw);
    plot([(rwtm-sbval)' ;(rwtm-sbval)'],[(muz-muzer); (muz+muzer)],'Color',ps.mucol,'Linewidth', ps.errlw);
    
    hold on;
    plot([0 max(rwtm-sbval)],[0,0],'k--','Linewidth',2)
    
    if ps.plotwn
         wnheight=1.1*abs(max(acz))+1.3;
        plotwnbnds(vls.wnon, vls.wnoff, vls.subvl, wnheight)
       
    
    end
    if (ps.plotruns)
       basht=1.1*abs(max(acz));
       shiftht=basht+.3;
       subht=basht+.6;
       initht=basht+.9
       revht=basht+1.2;
       
       %plot bs
       basruns=vls.basruns;
       bashtvc=basht*ones(length(basruns),1);
       xvec=rwtmall(basruns)-sbval;
       plot([xvec';xvec'],[bashtvc'-.1;bashtvc'+.1],'r');
       
       %plotshift
       shiftruns=vls.shiftruns;
       for ii=1:length(shiftruns)
           if(~isempty(shiftruns{ii}))
                shftvls=shiftruns{ii};
                shifthtvc=shiftht*ones(length(shftvls),1);
                xvec=rwtmall(shftvls)-sbval;
                plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'c');
           end
       end
           
       
       initruns=vls.initruns;
       for ii=1:length(initruns)
           if(~isempty(initruns{ii}))
                shftvls=initruns{ii};
                shifthtvc=initht*ones(length(shftvls),1);
                xvec=rwtmall(shftvls)-sbval;
                plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'k');
           end
       end

     subruns=vls.subruns;
       for ii=1:length(subruns)
           if(~isempty(subruns{ii}))
                shftvls=subruns{ii};
                shifthtvc=subht*ones(length(shftvls),1);
                xvec=rwtmall(shftvls)-sbval;
                plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'g');
           
       end
       end
      revruns=vls.revruns;
       for ii=1:length(revruns)
           if(~isempty(revruns{ii}))
                shftvls=revruns{ii};
                shifthtvc=revht*ones(length(shftvls),1);
                xvec=rwtmall(shftvls)-sbval;
                plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'m');
           end
       end
        
    end
   
    box off;
    title(bname )
    
    
    
 function []=plotwnbnds(wnon, wnoff, subvl, wnheight)
     
     for ii=1:length(wnon)
         for jj=1:length(wnon{ii})
            xvec1=wnon{ii}(jj)-subvl
            xvec2=wnoff{ii}(jj)-subvl
            wnht=wnheight*ones(length(xvec1),1);
            plot([xvec1;xvec1],[wnht'-.1; wnht'+.1],'k')
            plot([xvec2;xvec2],[wnht'-.1;wnht'+.1],'k')
            plot([xvec1;xvec2],[wnht';wnht'],'k')
         end
     end
          
          
          