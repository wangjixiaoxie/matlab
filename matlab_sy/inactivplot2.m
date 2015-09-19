figure
for nt=1:2
    for ii=1:length(timind)
        for jj=1:3
            clrvl='krb'
            linestyl{1}='-'
            linestyl{2}='-'
            linestyl{3}='--'
            lwid=[1 2 1];
            curtm=timind(ii);
            subplot(3,2,ii*2+nt-2)
            
            if(~isempty(mnhstout{nt}{curtm}{jj}))
                stairs(edges{nt}, mnhstout{nt}{curtm}{jj},'Color',clrvl(jj), 'Linestyle', linestyl{jj},'Linewidth',lwid(jj));
                hold on;    
                box off;
            
                nvl=num2str(nvls{nt}{curtm}{jj});
                cvl=num2str(cv{nt}{curtm}{jj});
                text(edges{nt}(jj*3),.2,[nvl ' ' cvl],'Color',clrvl(jj));
            end
        end
    end
end
