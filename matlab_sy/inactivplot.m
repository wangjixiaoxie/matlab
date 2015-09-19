figure
for ii=1:4
   
    subplot(4,2,2*ii-1)
    plotmeanstdev(avls.chkarray{1},avls.stimes,avls.etimes,graphvals,1,avls.chunkout{1})
    axis([13.2+ii 13.8+ii 2350 2800])
    subplot(4,2,2*ii)
    
    plotmeanstdev(avls.chkarray{2},avls.stimes,avls.etimes,graphvals,1,avls.chunkout{2})
    axis([13.2+ii,13.8+ii 6200 7000])
end


timind=[16 17]
edges{1}=[2350:35:2600]
edges{2}=[6200:70:7000]

for nt=1:2
    for ii=1:length(timind)
        curtm=timind(ii);
        tmon=graphvals.timon{curtm}
        tmoff=graphvals.timoff{curtm}
        dt=graphvals.date{curtm}
        
        
        
        tim1=datenum(dt,'yyyy-mm-dd');
        
        tim2=datenum([dt ' ' tmon],'yyyy-mm-dd HH:MM:SS')
        tim3=datenum([dt ' ' tmoff],'yyyy-mm-dd HH:MM:SS')
        tim4=ceil(tim3)
        tim2b=tim2+0.75/24;
        tim3b=tim3+0.75/24;
        
        vals=getvals(avls.fvcomb{nt},1,'TRIG');
        
        ind{1}=find(vals(:,1)>tim1&vals(:,1)<tim2);
        ind{2}=find(vals(:,1)>tim2b&vals(:,1)<tim3);
        ind{3}=find(vals(:,1)>tim3b&vals(:,1)<tim4);
        
        for kk=1:3
            if(~isempty(ind{kk}))
                chkvals{kk}=vals(ind{kk},2);
                hstout{nt}{curtm}{kk}=histc(chkvals{kk},edges{nt});
                mnhstout{nt}{curtm}{kk}=hstout{nt}{curtm}{kk}/sum(hstout{nt}{curtm}{kk})
                nvls{nt}{curtm}{kk}=length(chkvals{kk})
                stdev{nt}{curtm}{kk}=std(chkvals{kk})
                mnvls{nt}{curtm}{kk}=mean(chkvals{kk})
                cv{nt}{curtm}{kk}=std(chkvals{kk})/mean(chkvals{kk})
            else
                chkvals{kk}=[];
                
            end
        end
    end
end



%now to plot
figure
for nt=1:2
    for ii=1:length(timind)
        for jj=1:3
            clrvl='krk'
            linestyl{1}='-'
            linestyl{2}='--'
            linestyl{3}='-'
            curtm=timind(ii);
            subplot(3,2,ii*2+nt-2)
                
            stairs(edges{nt}, mnhstout{nt}{curtm}{jj},'Color',clrvl(jj), 'Linestyle', linestyl{jj},'Linewidth',2);
            hold on;    
            
            
            nvl=num2str(nvls{nt}{curtm}{jj});
            cvl=num2str(cv{nt}{curtm}{jj});
            
            text(edges{nt}(jj),.5,[nvl ' ' cvl],'Color',clrvl(jj));

        end
    end
end


        
        
           

