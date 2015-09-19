%rewritten 11.21.08

% plot_timecourse2.m

function [os] = plot_tmcourse2rev(sumbs, ps)

for ii=1:length(sumbs)
    ps=setplotvls(ps);
    ps.plotruns=0;
    ps.FILL=1
    ps.plotwn=0;
    if(isfield(ps,'adjx'))
        ps.adjx=ps.adjx;
    else
        ps.adjx=0;
    end
    ppf=ps.ppf
    plotcount=3;
    
    pn=mod(ii,ppf)
    %if this is the last plot
%      if(~exist('axin'))
%          if(pn==0)
%             pn=ppf;
%          end
%         if(pn==1)
%             figure;
%         end
%         ax(pn)=subplot(ppf,1,pn)
%      else
%          axes(axin);
%          ax=axin
%      end
    
    if(exist('plotind'))
        ps.plotind=plotind;
        [vls]=getvals(sumbs(ii),ps);
    else
        [vls]=getvals(sumbs(ii),ps)
    end
    if(isfield(ps,'STIM'))
        ps.mucol='r'
    else
        ps.mucol=[0.4 0.4 1]
    end
    [os]=plotvals(vls,ps);
end


    
    function [ps] = setplotvls(ps)
    
        ps.lw=1
        ps.msize=4
        ps.errlw=1
        
%         ps.mucol=[0.4 0.4 1]
%         if(exist('STIM'))
%             ps.mucol='r'
%         end
        
        ps.mufillcol=[0.82 0.82 0.82]
        ps.acfillcol=[0.92 0.96 0.98]
        ps.accol='k'
        ps.meanlw=2
        ps.errht=50
        ps.acmark='o'
        ps.mumark='o'
        ps.ppf=3;
  
  function   [vls]=getvals(sumbs,ps);
      vls.ntind=sumbs.ntind;
      ntind=sumbs.ntind;
      initind=ps.initind;
     
      if(isfield(sumbs,'rawtimes'))
          vls.subvl=(sumbs.rawtimes(initind,1));
      else
          if(isfield(sumbs,'tmvec'))
            vls.subvl=(sumbs.tmvec(initind,1));
          else
              vls.subvl=(sumbs.tmvc(initind,1));
              sumbs.tmvec=sumbs.tmvc;
          end
          end
      if(isfield(sumbs,'mulist'))
        indvl=find(sumbs.mulist~=0)
        if(isfield(ps,'plotind'))
            if(ps.plotind)
            indvl=ps.plotind
            end
        end
      
        for ii=1:length(indvl)
            ind=indvl(ii)
                
                    if(isfield(ps,'mulistind'))
                        if(ps.mulistind)
                            indvl2=sumbs.mulist(ind);
                        else
                            indvl2=ind
                        end
                    else
                        indvl2=ind
                    end
                    vls.acz(ii)=sumbs.acz(ntind,indvl2);
                    vls.acmean(ii)=sumbs.acmean(ntind,indvl2)/3000;
%                     vls.aczerr(ii)=sumbs.acerrz(ntind,indvl2)
%                     vls.acstdv(ii)=sumbs.acstdv(ntind,indvl2)/3000;
                    vls.muz(ii)=sumbs.muz(ntind,indvl2);
                    vls.mumean(ii)=sumbs.mumean(ntind,indvl2)'/3000
%                     vls.muzerr(ii)=sumbs.muerrz(ntind,indvl2)
%                     vls.mustdv(ii)=sumbs.mustdv(ntind,indvl2)/3000;
                
        end

    if(ps.mulistind)
        
        vls.rwtms=sumbs.tmvec(sumbs.mulist(indvl),1)
    else
        vls.rwtms=sumbs.tmvec(indvl,1)
    end
        vls.rawtmsall=sumbs.tmvec(:,1);
    vls.basruns=sumbs.basruns;
    vls.revruns=sumbs.revruns;
    vls.shiftruns=sumbs.shiftruns
    vls.initruns=sumbs.initind;
    vls.subruns=sumbs.subruns;
    vls.initmean=sumbs.initmean;
    vls.mulist=sumbs.mulist;
    vls.bname=sumbs.bname;
    vls.wnon=sumbs.flr_wnon
    vls.wnoff=sumbs.flr_wnoff
      else
        if(exist('plotind'))
            indvl=plotind
        else
            plotind=1:length(sumbs.tmvec(:,1))
        end

          vls.rwtms=sumbs.tmvec(plotind,1);
          [sortout,sortind]=sort(vls.rwtms);
          plotind=plotind(sortind);
          vls.rwtms=sumbs.tmvec(plotind,1);
          vls.rwtms=sumbs.tmvec(plotind);
          vls.muz=sumbs.muz(plotind);
          vls.acz=sumbs.acz(plotind);
          
          vls.acmean=sumbs.acmean(plotind)/3000;
          vls.mumean=sumbs.mumean(plotind)/3000
          vls.initmean=sumbs.initmean;
      end
  function [os]=plotvals(vls,ps);
    axes(ps.ax);
    rwtm=vls.rwtms
    acz=vls.acmean;
     muz=vls.mumean;
     
     if(ps.netdiff)
         os.netdiff=-(vls.acmean*1000-vls.mumean*1000);
         os.nettms=vls.rwtms-vls.subvl;
     end
     
     if(ps.plotdiff)
         acz=(vls.acmean*3000-vls.initmean)/3;
         muz=(vls.mumean*3000-vls.initmean)/3;
     end
     
     initmean=vls.initmean;
    if(isfield(vls,'rawtmsall'))
        rwtmall=vls.rawtmsall;
        sbval=vls.subvl;
        bname=vls.bname
%         aczer=vls.acstdv./sqrt(length(vls.acmean))
%         muzer=vls.mustdv./sqrt(length(vls.mumean));
    end
if(~exist('sbval'))
    sbval=rwtm(1);
    aczer=0
    muzer=0
end
  if(ps.FILL)
        xvec=rwtm-sbval
         if(ps.adjx)
            xvec=xvec-xvec(1)+1;
        end
        
        xvec=makerow(xvec)
        fillx=[xvec xvec(end:-1:1)]

        muvec=[acz]
        muvec=makerow(muvec);
        muvec2=[muz]
        muvec2=makerow(muvec2);
        filly=[muvec  muvec2(end:-1:1)]

% acvec=[rvacmn+rvacer]
% acvec2=[rvacmn-rvacer]
% filly2=[acvec acvec2(end:-1:1)]

        mufillcol=[0.82 0.82 0.82]  
% acfillcol=[0.92 0.96 0.98]

% %     mupts=mumean(1:length(indvl))
%     acpts=acmean(1:length(indvl))
%     yvec=[avls.initmean{notevl}*ones(length(mupts),1);mupts(end:-1:1)']
    

    hold on;
%     fill(xvec,yvec,acfillcol);
% filly([1 end])=0
% filly2([1 end])=0
    fill(fillx,filly,mufillcol,'edgecolor','w');
%     fill(fillx,filly2,mufillcol);
%     
    end  





    plot(xvec,acz,'Color',ps.accol,'Linewidth',ps.lw);
    hold on;
     plot(xvec,acz,'Linestyle','none','Marker',ps.acmark,'Color',ps.accol,'Markersize',ps.msize,'MarkerFaceColor',ps.accol); 
    plot(xvec,muz,'Color',ps.mucol,'Linewidth',ps.lw);
    plot(xvec,muz,'Linestyle','none','Marker',ps.mumark,'Color',ps.mucol,'Markersize',ps.msize,'MarkerFaceColor',ps.mucol);
    
    %plot vertical lines as error bars

%     plot([(xvec); (xvec)],[(acz-aczer); (acz+aczer)],'k','Linewidth',ps.errlw);
%     plot([(xvec) ;(xvec)],[(muz-muzer); (muz+muzer)],'Color',ps.mucol,'Linewidth', ps.errlw);
    
    hold on;
    plot([0 max(xvec)],[initmean(1)/3000,initmean(1)/3000],'k--','Linewidth',2)
    
    if ps.plotwn
         wnheight=1.1*abs(max(acz))+1.3;
        plotwnbnds(vls.wnon, vls.wnoff, vls.subvl, wnheight)
       
    
    end


    if(ps.plotarrow)
        
        for ii=1:length(xvec)
            origarrow=[xvec(ii) acz(ii) xvec(ii) muz(ii)];
            shiftdrxn=0;
            plotarrows(origarrow,ps);
        end
    end
        
%     
%     if (ps.plotruns)
%        basht=1.1*abs(max(acz));
%        shiftht=basht+.3;
%        subht=basht+.6;
%        initht=basht+.9
%        revht=basht+1.2;
%        
%        %plot bs
%        basruns=vls.basruns;
%        bashtvc=basht*ones(length(basruns),1);
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
        
        
                  
          