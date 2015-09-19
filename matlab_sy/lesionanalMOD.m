function [outvls]=lesionplot3(prevls,preind,pstind)
figure
figstoplot=[1 3]
daystoplot=[1 2 3  ]
postanaldy=3;
indtoplot=[1:5]
%   avlspath='/oriole5/evrenlesiondata/sumdata'
%   avlsmat='postshiftout.mat'
%     
%   cmd=['cd ' avlspath];
%   eval(cmd);
%   cmd=['load ' avlsmat];
%   eval(cmd);

load postshiftout.mat

edges=[2100:20:4800]
  %first plot hists. 
  axct=0;
 for indvl=1:length(indtoplot)
     crind=indtoplot(indvl);
     for dyvl=1:length(daystoplot)
         axct=axct+1;
         
         crdy=daystoplot(dyvl);
         ax(axct)=subplot(5,3,(indvl-1)*3+dyvl);
         crvls=prevls.crwnvls{crind}{crdy};
         hitind=find(crvls(:,3)==1);
         hstall=histc(crvls(:,2),edges)./length(crvls(:,2));
         hsthit=histc(crvls(hitind,2),edges)./length(crvls(:,2));
         if(dyvl==1)
             initmean=mean(crvls(:,2))/1000;
             
             hold on;
             inithsthit=hsthit;
         end
         stairs(edges/1000,hstall,'k')
         hold on;
         stairs(edges/1000,inithsthit,'Color',[0.8 0.8 0.8])
         stairs(edges/1000,hsthit,'r');
         plot([initmean initmean],[0 0.5],'k--')
     end
 end
 linkaxes(ax(1:3));
 linkaxes(ax(4:6));
 linkaxes(ax(7:9));
 linkaxes(ax(10:12));
         
         
         
         if(ismember(1,figstoplot))
             figure
    %ind should be 4.
    %exsong(1)
    %prelesion canonical
    exsong(1).ind=2
    exsong(1).TYPE='PRE'
%     exsong(1).ax=subplot(141);
    exsong(1).runind=1;
    exsong(1).day=3;
    exsong(1).daytoplot=1:4;
    %exsong(2)
    %postlesion canonical
    exsong(2).ind=2;
    exsong(2).TYPE='POS'
%     exsong(2).ax=subplot(142);
    exsong(2).runind=1;
    exsong(2).day=3;
    exsong(2).daytoplot=1:6;
%     exsong(3).ind=6;
%     exsong(3).TYPE='PRE'
%     exsong(3).ax=subplot(143);
%     exsong(3).runind=1;
%     exsong(3).day=3;
    %     %prelesion canonical
%     exsong(4).ind=10;
%     exsong(4).TYPE='PRE'
%     exsong(4).ax=subplot(144);
%     exsong(4).runind=1;
    
    %exsong(4)
    %postlesion canonical
%     exsong(4).ind=5;
%     exsong(4).TYPE='POS'
%     exsong(4).ax=subplot(248)
%     exsong(4).runind=1;
    precol='k'
    postcol{1}='r'
    postcol{2}='b'
    
    for ii=1:length(exsong)
          clear ps ax;
          ps.tst=1;
          
         birdind=exsong(ii).ind;
         runind=exsong(ii).runind;
%          ax(ii)=exsong(ii).ax;
         if(exsong(ii).TYPE=='PRE')
             ps.vlsprewn=prerawvls{birdind}.tst;
             ps.vlspostwn=prevls.crwnvls{birdind}{postanaldy}
             ps.allvls=prevls.crwnvls{birdind};
         else
            ps.vlsprewn=postrawvls{birdind}.tst{runind};
            ps.vlspostwn=postvls.crvls{birdind}{runind}{postanaldy}
            ps.allvls=postvls.crvls{birdind}{runind};
         end
          
%             [ps.THRESHVAL,ps.THRESHIND]=calc_escapethresh(ps.vlspostwn);
            %calc initescapes.
            
         
%             ps.prevls=vls{birdind}.tst;
%             ps.postvls=vls{birdind}.wn;
%             [ps.THRESHVAL,ps.THRESHIND]=calc_escapethresh(ps.postvls);
            %calcinitescapes
%             vlspstall=vlspstallout{birdind};
%             mindate=min(unique(floor(vlspstall(:,1))));
%             minind=find(floor(vlspstall(:,1))==mindate);
%             vlsinit=vlspstall(minind,:);
%             [ps.IN_THRESHVAL,ps.IN_THRESHIND]=calc_escapethresh(vlsinit);
%             ps.vlsinit=vlsinit;
%          end
         ps.edges=edges;
         ps.NORM=1;
         
         subplot(length(exsong),1,ii);
         ps.precol=precol;
         ps.postcol=postcol{ii};
         ps.daytoplot=exsong(ii).daytoplot;
         plotmnstdv(ps);
%         [outvls(ii)]=plotexhist(exsong(ii),ps);
%         axis([2.1 2.7 0 0.8])
    axis square;
    end

    
%     linkaxes(ax(1:length(exsong)));
    
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
  end
if(ismember(2,figstoplot))
    ps.STIM=1;
    ps.ax=subplot(244);
    for ii=1:length(postshiftout)
       
       
           postout(ii)=mean(postshiftout{ii})
      
    end
    
    combvls{1}{1}{1}.acshift=preshift.vals
    combvls{1}{1}{1}.mushift=postout;
    plotgroupbar(combvls,ps);
    axis square
    axis([2 5 -50 150])
end

if(ismember(3,figstoplot))
    figure
    preind=1:8
    pstind=[1:8 10 12:14]
    ps.ax=subplot(243);
    plot(preshift.cv_vals(preind),preshift.vals(preind),'ko','MarkerSize',4,'MarkerFaceColor','k')
    hold on;
    plot(postshift.cv_vals(pstind),postshift.vals(pstind),'ro','MarkerSize',4,'MarkerFaceColor','r')
    axis square
    axis([0 0.025 -20 150])
end

  function [outvls]=plotexhist(exsong,ps)
      axes(exsong.ax);
      hstpre=histc(ps.vlsprewn(:,2),ps.edges)./length(ps.vlsprewn);
      hstpst=histc(ps.vlspostwn(:,2),ps.edges)/length(ps.vlspostwn);
      
      hstesc=histc(ps.vlspostwn(ps.THRESHIND,2),ps.edges)/length(ps.vlspostwn);
%       if(isfield(ps,'vlsinit'))
%          hstinit=histc(ps.vlsinit(ps.IN_THRESHIND,2),ps.edges)/length(ps.prevls);
%           
%       end
      mnpre=mean(ps.vlsprewn(:,2))/1000;
      mnpst=mean(ps.vlspostwn(:,2))/1000;
      mnesc=mean(ps.vlspostwn(ps.THRESHIND),2)/1000;
      cvpre=std(ps.vlsprewn(:,2))./mean(ps.vlsprewn(:,2));
      cvpst=std(ps.vlspostwn(:,2))./mean(ps.vlspostwn(:,2));
      npre=length(ps.vlsprewn);
      npst=length(ps.vlspostwn);
      edgesout=ps.edges(1:end-1);
      hstpreout=hstpre(1:end-1);
      hstpstout=hstpst(1:end-1);
      hstesc=hstesc(1:end-1);
      credgediff=edgesout(2)-edgesout(1);
%       edgesout=edgesout+.5*credgediff;
      stairs(edgesout/1000,hstpreout,'k');
      hold on;
      stairs(edgesout/1000,hstpstout,'r');
      stairs(edgesout/1000,hstesc,'Color',[0.8 0.8 0.8]);
      if(isfield(ps,'vlsinit'))
          hstinit=hstinit(1:end-1);
         stairs(edgesout/1000,hstinit,'b');
      end
      
                    plot(mnpre,.4,'Marker','v','Color','k','MarkerFaceColor','k');
                    plot(mnpst,.4,'Marker','v','Color','r','MarkerFaceColor','r');
                    plot(mnesc,.4,'Marker','v','Color',[0.8 0.8 0.8]);
                    text(2.6,.7,['n=' num2str(npre)],'Color','k');
                    text(2.6,.65,['mean=' num2str(mnpre)],'Color','k');
                    text(2.6,.6,['cv=' num2str(cvpre)],'Color','k');
                    text(2.6,.55,['n=' num2str(npst)],'Color','r');
                    text(2.6,.5,['mean=' num2str(mnpst)],'Color','r');
                    text(2.6,.45,['cv=' num2str(cvpst)],'Color','r');
                    outvls.mnpst=mnpst;
                    outvls.cvpst=cvpst;
                    outvls.mnpre=mnpre;
                    outvls.cvpre=cvpre;
                    [h,p]=ttest2(ps.vlsprewn(:,2),ps.vlspostwn(:,2))
                    outvls.pval=p;
      
      
      %text describing mnvl cv
      function []= plotmnstdv(ps)
          unqdays=unique(floor(ps.vlsprewn(:,1)))
         
          for ii=1:length(unqdays)
              crind=find(floor(ps.vlsprewn(:,1))==unqdays(ii));
              crmn(ii)=mean(ps.vlsprewn(crind,2));
              crstd(ii)=std(ps.vlsprewn(crind,2));
              crx(ii)=ii-length(unqdays)
          end
          initmn=mean(crmn);
           for ii=1:length(unqdays)
            plot(crx(ii),initmn-crmn(ii),'Marker','o','MarkerEdgeColor',ps.precol,'MarkerFaceColor',ps.precol);
              hold on;
              plot([crx(ii) crx(ii)],[initmn-crmn(ii)+crstd(ii) initmn-crmn(ii)-crstd(ii)],ps.precol);
          end
          
          clear crx crmn crstd
          for ii=1:length(ps.daytoplot)
              crvls=ps.allvls{ii};
              crmn=mean(crvls(:,2));
              crmnout(ii)=crmn;
              crstd=std(crvls(:,2));
              plot(ii,initmn-crmn,'Marker','o','MarkerEdgeColor',ps.postcol,'MarkerFaceColor',ps.postcol);
              hold on;
              plot([ii ii],[initmn-crmn+crstd initmn-crmn-crstd],ps.postcol);
          end
          tst=1;