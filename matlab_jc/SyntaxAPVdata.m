% SyntaxAPVdata

% 1. Label files and create a batch of *.cbin files (e.g. 'batchfiles')
% 2. Create avls structure for ACSF runs
ptsynanal_r39g39.m




% 3. Analyze and plot data
        % Calculate learning curves
            avls1=avls;
            avls1.pvls=avls.pvls(avls.indCurve);
            avls1.cvl=avls.cvl(avls.indCurve);
            avls1.MKSYN=1:length(avls1.pvls);
            syn_settimestampslabelsCURVE(avls1) % saved as /datsum/sumdataCURVES.mat
        % Calculate learning points + CIs
            avls2=avls;
            avls2.pvls=avls.pvls(avls.indPoint);
            avls2.cvl=avls.cvl(avls.indPoint);
            avls2.MKSYN=1:length(avls2.pvls);    
            syn_settimestampslabelsPOINT(avls2) % saved as /datsum/sumdataPOINTS.mat
        % Plot learning curves    
        load sumdataCURVES.mat
        clear times probs probs05 probs95
        for ii=1:length(tms)
            times{avls.indCurve(ii)}=tms{ii};
            probs{avls.indCurve(ii)}=pout{ii};
            probs05{avls.indCurve(ii)}=p05out{ii};
            probs95{avls.indCurve(ii)}=p95out{ii};
        end
        load sumdataPOINTS.mat
        for ii=1:length(tms)
            times{avls.indPoint(ii)}=tms{ii};
            probs{avls.indPoint(ii)}=pout{ii};
            probs05{avls.indPoint(ii)}=p05out{ii};
            probs95{avls.indPoint(ii)}=p95out{ii};
        end

        figure;hold on;
            for ii=avls.indCurve
                unq_days=unique(floor(times{ii}));
                for jj=1:length(unq_days)
                crday=unq_days(jj);
                crdayind=find(floor(times{ii})==crday);
                xvls=[times{ii}(crdayind),times{ii}(crdayind(end:-1:1))];
                yvls=[probs05{ii}(crdayind) probs95{ii}(crdayind(end:-1:1))];
                fill(xvls,yvls,'b')
                end
            end   
        % Plot learning points + CIs
        for ii=avls.indPoint
                xvls=[times{ii};times{ii}(end:-1:1)];
                yvls=[probs05{ii} probs95{ii}(end:-1:1)]; 
                if isempty(find(avls.indPointAC==ii))
                    plot(xvls,yvls,'r-','Linewidth',2)
                    plot(mean(xvls),mean(yvls),'r.','Markersize',15)
                else
                    plot(xvls,yvls,'b-','Linewidth',2)
                    plot(mean(xvls),mean(yvls),'b.','Markersize',15)
                end
        end
        t1=times{1};
        if length(times{end})==1;t2=times{end};else t2=max(times{end});end
        for i=1:length(avls.indPointAC); pp(i)=probs{avls.indPointAC(i)}; end
        if length(pp)==1;p1=pp;else p1=mean(pp);end
        plot([t1 t2],[p1 p1],'k')

