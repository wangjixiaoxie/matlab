%rewritten 12.10.09
%so as to plot acsf and muscimol including all aclist runs.
%(even days when no muscimol)

% plot_timecourse2.m

%THis version explicitly includes beginning and end of day in different
%colors

function [netdiff,vls] = plot_tmcourse5MOD2(crbs, ps,avls,vls)
    

    ps=setplotvls(ps);
    ps.plotruns=0;
%     ps.FILL=1
    ps.plotwn=0;
    if(isfield(ps,'adjx'))
        ps.adjx=ps.adjx;
    else
        ps.adjx=0;
    end
    ppf=ps.ppf
    plotcount=3;
    
    
  
    [vls]=getvals(crbs,ps,avls);
   
    if(isfield(ps,'STIM'))
        ps.mucol='r'
    else
        ps.mucol=[0.4 0.4 1]
    end
%     ps.mucol=[0.4 0.4 1]
ps.mucol='r'
    [netdiff]=plotvals(vls,ps,avls);


    
    function [ps] = setplotvls(ps)
    
        ps.lw=1
        ps.ntind=1;
        ps.msize=4
        ps.errlw=1
        
%         ps.mucol=[0.4 0.4 1]
%         if(exist('STIM'))
%             ps.mucol='r'
%         end
        
ps.mufillcol=[0.92 0.96 0.98]
 ps.acfillcol=[.98 .8 .9]
%         ps.acfillcol=[0.92 0.96 0.98]
        ps.accol='k'
        ps.meanlw=2
        ps.errht=50
        ps.acmark='o'
        ps.mumark='o'
        ps.ppf=3;
  
  %this version is rewritten to calculate beginning and endvls
        function   [vls]=getvals(crbs,ps,avls)
      vls.ntind=ps.ntind;
      ntind=ps.ntind;
      acindvl=ps.acindvl
      muindvl=find(crbs.mulist(ps.muindvl)~=0);
      muindvl=ps.muindvl(muindvl);
      if(ps.STIM)
          crbs.tmvc=crbs.tmvec
      end
      vls.zerotime=floor(crbs.tmvc(ps.startind,1))+.5;
        for ii=1:length(acindvl)
                crind=acindvl(ii)
        end
        
        for ii=1:length(muindvl)
            
            crmuind=muindvl(ii)
            cr_murun=crbs.mulist(crmuind);
            vls.rawmutms{ii}=avls.adjvls{1}{cr_murun}(:,1)-vls.zerotime;  
            vls.rawmuvls{ii}=avls.adjvls{1}{cr_murun}(:,2)/3000;
            vls.mumean(ii)=crbs.mumean(1,cr_murun)./ps.divfac;
            vls.mudiff(ii)=1000*((vls.mumean(ii)-crbs.initmean(1)./ps.divfac));
            vls.mustdv(ii)=crbs.mustdv(1,cr_murun)./ps.divfac;
            
            vls.mutms(ii)=mean(crbs.tmvc(cr_murun,:));
                 
        end
        %this is currently written only to use second ac value, need both!
        for ii=1:length(ps.runstoplot)
            crind=ps.runstoplot(ii);
           indvl2=crbs.mulist(ii);
                        
                        %ACINDTOPLOT is 1 or 2, or both depending on which
                        %AC columns to include
                        for indnum=ps.ACINDTOPLOT   
                            acind(indnum)=crbs.aclist(crind,indnum);
                            [tmtmp(indnum),mntmp(indnum),stdtmp(indnum)]=mk_acmean(avls,acind(indnum))
                            
                        end
                        calcmn=nanmean(mntmp);
                        calcstd=nanmean(stdtmp);
                        calctm=nanmean(tmtmp);
                            vls.acmean(ii)=calcmn/ps.divfac;
                            vls.acdiff(ii)=1000*((calcmn-crbs.initmean(1))./ps.divfac);
                            vls.acstdv(ii)=calcstd./ps.divfac;
                            vls.actms(ii)=calctm;
                            vls.initmean=crbs.initmean(1);
                            vls.initsd=crbs.initsd(1);
                            vls.acz(ii)=(calcmn-crbs.initmean(1))./crbs.initsd(1)
       end
           
        
    
    vls.mutms=vls.mutms-vls.zerotime;
%     vls.init_tms=vls.initrwtms-vls.zerotime;
%     vls.end_tms=vls.endrwtms-vls.zerotime;
    vls.basmean=crbs.initmean./ps.divfac;
%   
            function [out_tm,out_mn,out_std]=mk_acmean(avls,acind)
                outptvls=[];
                outtmvls=[];
                acindvl=find(acind>0);
                acind=acind(acindvl);
                
                for ii=1:length(acind)
                    outptvls=[outptvls ;avls.adjvls{1}{acind(ii)}(:,2)];
                    outtmvls=[outtmvls ;avls.adjvls{1}{acind(ii)}(:,1)];
                end
                out_tm=mean(outtmvls);
                out_mn=mean(outptvls);
                out_std=std(outptvls);    
                
                    
            function [os]=plotvals(vls,ps,avls)
  %this function has to plot muscimol vals
  %this function has to plot muscimol line
  %this function has to plot ac vls (two colors)
  %this function has to plot acline
  %this function has to plot two fills.
%                 
%                 axes(ps.ax);
%     rwtm=floor(vls.rwtms)
    acz=vls.acmean;
     muz=vls.mumean
     nzind=find(vls.mumean)

hold on;
if(ps.plotz)
    combacvls=vls.acz;
    muvls=(vls.mumean-vls.initmean)./vls.initsd;
    stvls=vls.mustdv;
elseif(ps.plotdiff)
    combacvls=vls.acdiff;
    muvls=vls.mudiff;
    stvls=1000*vls.mustdv
    
else
   combacvls=vls.acmean; 
    muvls=vls.mumean;
end
if(ps.flip)
    combacvls=-combacvls;
    muvls=-muvls;
    vls.acz=-vls.acz
    vls.acdiff=-vls.acdiff;
    vls.mudiff=-vls.mudiff;
end


if(isfield(ps,'muindvl'))
    ps.mulistind=ps.muindvl;
    ps.muind=ps.muindvl;
end
    combactms=vls.actms(1:length(ps.mulistind))-vls.zerotime;
    
%     initactms=vls.init_tms;
%     
%     endactms=vls.end_tms;
    
    mutms=vls.mutms;
    
        xvec=combactms;
        xvecmu=makerow(mutms);
        xvec=makerow(xvec);

%   if(ps.FILL)
        
        fillxvecmu=[combactms(1) xvecmu combactms(end)]
        xvecbas=[min(fillxvecmu) max(fillxvecmu)]
        fillx=[xvec fillxvecmu(end:-1:1)]
        fillx2=[xvecmu  xvecmu(end) xvecmu(1)]

        muvec=[combacvls(1:length(ps.mulistind))]
        muvec=makerow(muvec);
        %for fill against acsf
        
        muvectop=makerow([muvls(1) muvls muvls(end)])
        muvecbot=makerow(muvls);
      
        muvecbas=[0 0]
                
        
        filly=[muvec  muvectop(end:-1:1)]
        filly2=[muvecbot muvecbot(1) muvecbot(1) ];


        mufillcol=ps.mufillcol;
acfillcol=ps.acfillcol;
    
if(ps.FILL)
    hold on;
    fill(fillx,filly,mufillcol,'edgecolor','none');
    fill(fillx2,filly2,acfillcol,'edgecolor','none');
%     
  end  
  start_time=[floor(avls.rawtimes(ps.startind,1))+.5] ;

  if isfield(ps,'plotwn')
      if(ps.plotwn)
   plotwn(avls,start_time,ps); 
  
      end
  end

  
  
% axis([6 15 -1 5])
if(isfield(ps,'plotlman'))
   if(ps.plotlman)
       xvecmu=xvec(nzind);
       xvecac=xvec
      [diffx,diffy]=getdiff(xvecmu,xvecac,muz(nzind),acz)
      hold on;
      plot(diffx(nzind),diffy(nzind),'Marker',ps.acmark','Color',ps.accol,'MarkerFaceColor',ps.accol,'Linewidth',ps.lw)
   end
end

if(ps.plotsep)
    
    if(ps.plotinitendline)
            for ii=1:length(initacvls)
            plot([initactms(ii) endactms(ii)],[initacvls(ii) endacvls(ii)],'Color',ps.accol,'Linewidth',ps.lw);
            hold on;
        end
        plot(initactms,initacvls,'Linestyle','none','Marker',ps.acmark,'Color',ps.accol,'Markersize',ps.msize,'MarkerFaceColor',ps.accol);
      plot(endactms,endacvls,'Linestyle','none','Marker',ps.acmark,'Color',ps.accol,'Markersize',ps.msize,'MarkerFaceColor',ps.accol); 
    end
      if(ps.plotmumarker)
        plot(xvecmu,ps.muht*ones(1,length(xvecmu)),'Marker','v','Color',ps.mucol,'MarkerFaceColor',ps.mucol,'Linestyle','none');  
      end
      if(ps.plotmuline)
          
          %this is in case mu run is not on last day before wn cameon
          if(isfield(ps,'plotmudash'))
              
             indlo=find(xvecmu<ps.plotmudashend);
             indhi=find(xvecmu>ps.plotmudashend);
             [initxvec]=[xvecmu(indlo) ps.plotmudashend];
             [initmu]=[muvls(indlo) muvls(indlo(end))];
             [finxvec]=[ps.plotmudashend xvecmu(indhi)] 
             [finmu]=[muvls(indlo(end)) muvls(indhi)]
            
%              plot(xvecmu,muvls,'Color',ps.mucol,'Linewidth',ps.lw);
             plot(initxvec, initmu,'Color',ps.mucol,'Linewidth',ps.lw,'Linestyle','--');
             plot(finxvec, finmu,'Color',ps.mucol,'Linewidth',ps.lw);
             plot(xvecmu,muvecbot,'Linestyle','none','Marker',ps.mumark,'Color',ps.mucol,'Markersize',ps.msize,'MarkerFaceColor',ps.mucol);
          else
            plot(xvecmu,muvls,'Color',ps.mucol,'Linewidth',ps.lw);
            plot(xvecmu,muvecbot,'Linestyle','none','Marker',ps.mumark,'Color',ps.mucol,'Markersize',ps.msize,'MarkerFaceColor',ps.mucol);
               plot([xvecmu ;xvecmu], [muvls+stvls;muvls-stvls],'Color',ps.mucol,'Linewidth',ps.lw)
          
          end
      end
      if(ps.plotrawmu)
         for ii=1:length(vls.rawmuvls)
             plot(vls.rawmutms{ii},vls.rawmuvls{ii},'o','MarkerFaceColor',ps.mucol,'MarkerEdgeColor',ps.mucol,'MarkerSize',ps.mumarksize);
         end
      end
    if(ps.plotallac)
        crvlcomb=[];
        crtmcomb=[];
        if (isfield(ps,'runstoplot'))
            nozeroind=1:length(ps.runstoplot)
        else
            nozeroind=find(avls.aclist(:,1)~=0);
        end
        for ii=1:length(nozeroind)
            crind=ii;
            if(ps.plotz)
                crvl=vls.acz(crind);
            elseif(ps.plotdiff)
                crvl=vls.acdiff(crind);
                crvlstd=1000*vls.acstdv(crind);
            else
                crvl=vls.acmean(crind);
            end
                crtm=vls.actms(crind);
                crvlcomb=[crvlcomb crvl]
%                 crvlsstd=vls.acstdv(cr
%                 crvlstd=vls.acstdv(crind);
                crtmcomb=[crtmcomb crtm];
            hold on;  
            if(ps.plotacline)
            plot(round(crtm-start_time),crvl,'o','MarkerFaceColor',ps.accol,'MarkerEdgeColor',ps.accol,'MarkerSize',ps.msize)
            plot([round(crtm-start_time) round(crtm-start_time)],[crvl-crvlstd  crvl+crvlstd],'k','Linewidth',ps.lw);
            plot(round(crtmcomb-start_time),crvlcomb,'k','Linewidth',ps.lw);
            end
if(ps.plotallacpts)
            if(avls.aclist(crind,1))
                 crvls=avls.allvls{1}{avls.aclist(crind,1)}
                hold on;       
                 plot(crvls(:,1)-start_time,crvls(:,2)/3000,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1.5)
                
            end
            if(avls.aclist(crind,2))
                 crvls=avls.allvls{1}{avls.aclist(crind,2)}
                hold on;       
                 plot(crvls(:,1)-start_time,crvls(:,2)/3000,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1.5)
                
            end
end
        end
        
        end
    

      if(ps.plotrawac)
        for ii=1:length(vls.initrwtms)
            plot(vls.initrwtmsall{ii},vls.crinitvls{ii},'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1.5);
            plot(vls.endrwtmsall{ii},vls.crendvls{ii},'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1.5);
        end
    end

else

    plot(xvec,combacvls,'Color',ps.accol,'Linewidth',ps.lw);
    hold on;
     plot(initactms,initacvls,'Linestyle','none','Marker',ps.acmark,'Color',ps.accol,'Markersize',ps.msize,'MarkerFaceColor',ps.accol);
      plot(endactms,endacvls,'Linestyle','none','Marker',ps.acmark,'Color',ps.accol,'Markersize',ps.msize,'MarkerFaceColor',ps.accol); 
    plot(xvecmu,muvls,'Color',ps.mucol,'Linewidth',ps.lw);
    plot([xvecmu xvecmu], [stdmuvls stdmuvls],'Color',ps.mucol,'Linewidth',ps.lw)
    plot(xvecmu,muvecbot,'Linestyle','none','Marker',ps.mumark,'Color',ps.mucol,'Markersize',ps.msize,'MarkerFaceColor',ps.mucol);
end   
 
    hold on;
%     plot([0 max(endactms)],[vls.basmean,vls.basmean],'k--','Linewidth',2)
    
    if ps.plotwn
         wnheight=1.1*abs(max(acz))+1.3;
        plotwnbnds(vls.wnon, vls.wnoff, vls.subvl, wnheight)
       
    
    end
    if(ps.netdiff)
    tmpind=find(avls.mulist(ps.runstoplot(nozeroind))~=0);
    matchind=nozeroind(tmpind);
    os.netdiff=(vls.acdiff(matchind)-vls.mudiff);
    os.mudiff=vls.mudiff;
    os.nettms=vls.actms(matchind)-start_time;
    
    end
    if(isfield(ps,'wnind'))
        axes(ps.wnax);
        ind=find(ps.wnind>0)
        %make xvec
        %make yvec
        for ii=1:length(ind)
           xvecout(2*ii-1)=xvec(ind(ii))-.5
           xvecout(2*ii)=xvec(ind(ii))+.5
           wnout(2*ii-1)=ps.wnbnds(ii);
           wnout(2*ii)=ps.wnbnds(ii);
        end
        xcomb=[xvecout xvecout(end:-1:1)]
        yout=[wnout ps.wnbot*ones(1,length(wnout))];
        fill(xcomb,yout,'r','EdgeColor','none')
    end


     function []=plotwn(avls,start_time,ps)
        ptvls=avls.wnptvls{1};
        tms=avls.wntms{1}-start_time;
        tmsout=[];
        ptout=[];
        for ii=1:length(tms)
            if(ii>1)
                tmsout=[tmsout tms(ii) tms(ii)];
                ptout=[ptout ptvls(ii-1) ptvls(ii)];
            else
                tmsout=[tmsout tms(ii) ];
                ptout=[ptout ptvls(ii)]
            end
        end
        
        xvls1=[makerow(tmsout) 9]; ;
        xvls2=xvls1(end:-1:1);
        yvls1=[makerow(ptout) ptout(end)]
        yvls1=yvls1/3000;
        yvls2=ps.wnfloor*(ones(length(xvls1),1));
%         
        xcomb=[makerow(xvls1) makerow(xvls2)]
        ycomb=[makerow(yvls1) makerow(yvls2)];
        if(ps.plotwnline==0) 
            fill(xcomb,ycomb,ps.wncol,'edgecolor','none');
        else
            plot(xvls1,yvls1,'Linestyle','--','Color',ps.wncol);
        end
        
    



   function [xdiff,ydiff]=getdiff(x1, x2, y1, y2)
       minvl=min([x1 x2])
       maxvl=max([x1 x2])
       tms=minvl:1:maxvl
       interpy1=interp1(x1,y1,tms);
       interpy2=interp1(x2,y2,tms);
       xdiff=tms;
       ydiff=interpy2-interpy1;
     
    
    
    
%        xvec=rwtmall(basruns)-sbval;
%        plot([xvec';xvec'],[bashtvc'-.1;bashtvc'+.1],'r');
%        
%        %plotshift
%        shiftruns=vls.shiftruns;
%        for ii=1:length(shiftruns)
%            if(~isempty(shiftruns{ii}))
%                 shftvls=shiftruns{ii};
%                 shifthtvc=shiftht*ones(length(shftvls),1);
%                 xvec=rwtmall(shftvls)-sbval;
%                 plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'c');
%            end
%        end
%            
%        
%        initruns=vls.initruns;
%        for ii=1:length(initruns)
%            if(~isempty(initruns{ii}))
%                 shftvls=initruns{ii};
%                 shifthtvc=initht*ones(length(shftvls),1);
%                 xvec=rwtmall(shftvls)-sbval;
%                 plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'k');
%            end
%        end
% 
%      subruns=vls.subruns;
%        for ii=1:length(subruns)
%            if(~isempty(subruns{ii}))
%                 shftvls=subruns{ii};
%                 shifthtvc=subht*ones(length(shftvls),1);
%                 xvec=rwtmall(shftvls)-sbval;
%                 plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'g');
%            
%        end
%        end
%       revruns=vls.revruns;
%        for ii=1:length(revruns)
%            if(~isempty(revruns{ii}))
%                 shftvls=revruns{ii};
%                 shifthtvc=revht*ones(length(shftvls),1);
%                 xvec=rwtmall(shftvls)-sbval;
%                 plot([xvec';xvec'],[shifthtvc'-.1;shifthtvc'+.1],'m');
%            end
%        end
% %         
%     end
   
    box off;
%     title(bname )
  
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
          

    function []=plotarrows(origarrow,ps)
% %         axes(axhist);
% %         if(shiftdrxn==0)
% %             arrow_col=[1 .6 .6]
% %         else
% %             arrow_col=[.6 .6 1]
% %         end
if(isfield(ps,'arrow_col'))
    arrow_col=ps.arrow_col;
    
else
    if(origarrow(4)-origarrow(2)<0)
        arrow_col=[1 .6 .6]
    else
        arrow_col=[.6 .6 1]
    end
    
end
        origlen=origarrow(4)-origarrow(2);
        if(origlen<0)
            y1=origarrow(2);
            
        else
            y1=origarrow(2);
            
        end
        y2=origarrow(4);
%         shiftlen=shiftarrow(3)-shiftarrow(1);
%         headwidth(1)=.8/abs(origlen);
% %         headwidth(2)=.01/abs(shiftlen);
%         headht(1)=.01/abs(origlen)
%         headht(2)=.01/abs(shiftlen)
axis square
hold on;
if (isfield(ps,'axbnds'))
    axis(ps.axbnds)
end
    
%             global ColorOrder, ColorOrder=[];
%             set(gca,'ColorOrder',arrow_col)
%                global LineWidthOrder, LineWidthOrder=[2];
            
            h=arrow([origarrow(1) y1],[origarrow(3) y2],'Length',3);
            set(h,'FaceColor',arrow_col);
            set(h,'EdgeColor','none');
            set(h,'Linewidth',2)
%         arrowh([origarrow(1) origarrow(3)],[y1 y2],arrow_col,[100 100]);
%         plot([origarrow(1) origarrow(3)],[origarrow(2) origarrow(4)],'Color',arrow_col,'Linewidth',2)
        
%         plot_arrow(shiftarrow(1),shiftarrow(2),shiftarrow(3),shiftarrow(4),'Color',arrow_col,'Linewidth',3,'headwidth',headwidth(2),'headheight',headht(2),'facecolor',arrow_col,'edgecolor',arrow_col)
        
%         axis equal;
        
 function []=plotedge(ybnds,xvecin)
%            minx=min([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            maxx=max([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            st_x=floor(minx);
%            max_x=ceil(maxx);
           
           for ii=1:length(xvecin)
               crx=xvecin(ii);
              xbnds(1)=crx-9/24;
              xbnds(2)=crx;
              xvec=[xbnds xbnds(end:-1:1)]
              yvec=[ybnds(1) ybnds(1) ybnds(2) ybnds(2)]
              fill(xvec,yvec,[0.7 0.7 0.7],'edgecolor','none');
              hold on;
               
               
           end       
                  
          