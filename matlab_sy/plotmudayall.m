%plotmuscdayscript
%a hack to write plotmeanstdev for a single day.

nnotes=3
% btlst{1}=[2 6 14 19 23]
 btlst{1}=[ 26 29 32 35 39 42]
btlst{2}=btlst{1}+1;
btlst{3}=btlst{1}+2;
refvals=[7094 2200 1552]
wnon=graphvals.wnon
    subval=floor(avls.chkarray{1}{1}(1,1));
    stimes=avls.stimes;
    etimes=avls.etimes;
%put in the white noise
    fillind=[1 .9 .9]
    mnvl=avls.mnvl;
    stdv=avls.stdv;


    
for ntcnt=1:3
for ii=1:length(btlst{1})
   
       
    subplot(length(btlst{1}),3,(3*(ii)+(ntcnt-3)))
    if(ntcnt==1)
        title(avls.pvls{btlst{2}(ii)})
        
%         for kk=1:length(wnon)
%             x1=datenum([wnon{kk}],'yyyy-mm-dd  HH:MM:SS')-subval;
%             x2=datenum([wnoff{kk}],'yyyy-mm-dd  HH:MM:SS')-subval;
%             y1=minbnds{kk};
%             y2=maxbnds{kk};
%             fill([x1 x2 x2 x1], [y1 y1 y2 y2], fillind,'EdgeColor','none')
%            
%            hold on; 
%         
%         
%     end
    
    
    end
        hold on;
    for tbinvl=1:3
         btvl=btlst{tbinvl}(ii);
        
    
    

   
          x=mean([stimes ;etimes],1)
          if tbinvl==1
              x=x1+0.5/24;
          elseif tbinvl==2
              x=x1+1/24;
          elseif tbinvl==3
              x=x1+1.5/24;
          end
              y=mnvl{ntcnt}(btvl);
            e=stdv{ntcnt}(btvl);
          txtht=graphvals.edges{ntcnt}(end)
          
          xvec=[x;x]
              yvec=[y+e;y-e]
              yvec2=[y;y]
             
%                 ofst=graphvals.offset(ii)
             
              
             if(tbinvl==2)
                 colvl=graphvals.col(1);
             else
                colvl=graphvals.col(2);
             end

                 plot(xvec, yvec, 'Color',colvl,'Linewidth',1);
              hold on;
             plot(xvec(1,:), yvec2(1,:), 'o','Color',colvl,'MarkerSize',2)
             box off;
             text(x, txtht, num2str(avls.nnote{ntcnt}{btvl}),'Fontsize',10)
                          plot([x1+0.5/24 x1+1.5/24], [refvals(ntcnt) refvals(ntcnt)],'b--','Linewidth',1)
             ylim([graphvals.edges{ntcnt}(1) graphvals.edges{ntcnt}(end)])
             xlim([x1+0.25/24 x1+1.75/24])  
            axis off;
               
             
          end
       end
end
    
    
