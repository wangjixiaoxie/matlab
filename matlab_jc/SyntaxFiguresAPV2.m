%%%% r37g7 - syntax
load /cardinal4/SyntaxAPV/r37summary.mat
    % Learning
        figure;hold on;
        % Experiment 1 - hit A, starting at baseline   %
        % runanalysis(1-(LClearnA1=='a'),1,1-LCbaseA,0.005,0)
            unidays=find(diff(TClearnA1)>8);
            unistart=[1 unidays+1];
            uniend=[unidays length(TClearnA1)];
            for i=1:length(unistart)
               unqtms=unique(TClearnA1(unistart(i):uniend(i)));
               times1=[];pmid=[];p05=[];p95=[];
               for tmind=1:length(unqtms)
                   selind=find(TClearnA1==unqtms(tmind));
                   times1(tmind)=unqtms(tmind);
                   pmid(tmind)=mean(e1.pmid(selind));
                   p05(tmind)=mean(e1.p05(selind));
                   p95(tmind)=mean(e1.p95(selind));
               end
                xvls=[times1(1:end),times1(end:-1:1)];
                yvls=[1-p05(1:end),1-p95(end:-1:1)];
                fill(xvls,yvls,'b')
                plot(times1,1-pmid,'r')
            end
            xv=[mean(timing3(Experiment1(5).fvAcatch)) mean(timing3(Experiment1(7).fvAcatch)) mean(timing3(Experiment1(4).fvAcatch)) 0 mean(timing3(Experiment1(9).fvAcatch))];
            xv=[2470 2480 2490 xv(4:5)];
            for i=[2 5]
                xvls=[xv(i) xv(i)];
                yvls=[probs05{i} probs95{i}];
                plot(xvls,yvls,'r-','Linewidth',2)
                plot(mean(xvls),mean(yvls),'r.','Markersize',15)
            end
            for i=[1 3]
                xvls=[xv(i) xv(i)];
                yvls=[probs05{i} probs95{i}];
                plot(xvls,yvls,'b-','Linewidth',2)
                plot(mean(xvls),mean(yvls),'b.','Markersize',15)
            end
            plot([min(TClearnA1)-50 max(TClearnA1)+20],[mean([probs{1} probs{3}]) mean([probs{1} probs{3}]) ],'k')
            ylim([0 0.5])
            xlim([2460 2555])
            
            
%%%% pu67bk2 - syntax
load /cardinal4/SyntaxAPV/pu67summary.mat
    % Learning
    figure;hold on;
        % Experiment 1 - hit A, starting at baseline   %
        % runanalysis(1-(LClearnA1=='a'),1,1-LCbaseA,0.005,0)
            unidays=find(diff(TClearn)>8);
            unistart=[1 unidays+1];
            uniend=[unidays length(TClearn)];
            for i=1:length(unistart)
               unqtms=unique(TClearn(unistart(i):uniend(i)));
               times1=[];pmid=[];p05=[];p95=[];
               for tmind=1:length(unqtms)
                   selind=find(TClearn==unqtms(tmind));
                   times1(tmind)=unqtms(tmind);
                   pmid(tmind)=mean(e1.pmid(selind));
                   p05(tmind)=mean(e1.p05(selind));
                   p95(tmind)=mean(e1.p95(selind));
               end
                xvls=[times1(1:end),times1(end:-1:1)];
                yvls=[1-p05(1:end),1-p95(end:-1:1)];
                fill(xvls,yvls,'b')
                plot(times1,1-pmid,'r')
            end
            res_vec=(LCcatch1=='a');
            guess=[];
                for i=1:1000
                    randsam=ceil(length(res_vec)*rand(1,length(res_vec)));
                    guess(i)=mean(res_vec(randsam));
                end
            xvls=([mean(timing3(Experiment1(1).fvAcatch)) mean(timing3(Experiment1(1).fvAcatch))]);
            yvls=([prctile(guess,5) prctile(guess,95)]);  
            plot(xvls,yvls,'b-','Linewidth',2)
            plot(xvls(1),mean(res_vec),'b.','Markersize',15)
            res_vec=(LCcatch2=='a');
            guess=[];
                for i=1:1000
                    randsam=ceil(length(res_vec)*rand(1,length(res_vec)));
                    guess(i)=mean(res_vec(randsam));
                end
            xvls=([mean(timing3(Experiment1(2).fvA)) mean(timing3(Experiment1(2).fvA))]);
            yvls=([prctile(guess,5) prctile(guess,95)]);  
            plot(xvls,yvls,'r-','Linewidth',2)
            plot(xvls(1),mean(res_vec),'r.','Markersize',15)
            for k=[2 4 6 8]
                res_vec=(LClearn_apv(k).data=='a');
                guess=[];
                for i=1:1000
                    randsam=ceil(length(res_vec)*rand(1,length(res_vec)));
                    guess(i)=mean(res_vec(randsam));
                end
                xvls=([mean(timing3(Experiment1(k).fvAcatch)) mean(timing3(Experiment1(k).fvAcatch))]);
                yvls=([prctile(guess,5) prctile(guess,95)]);
                plot(xvls,yvls,'r-','Linewidth',2)
                plot(xvls(1),mean(res_vec),'r.','Markersize',15)
            end
            plot([4235 4345],[mean(LCcatch1=='a') mean(LCcatch1=='a')],'k')
            
            
%%%% r39g39 - syntax
      % edit SyntaxAPVdata.m
        load /cardinal3/r39g9/datsum/avls.mat
        load /cardinal3/r39g9/datsum/sumdataCURVES.mat
        clear times probs probs05 probs95
        for ii=1:length(tms)
            times{avls.indCurve(ii)}=tms{ii};
            probs{avls.indCurve(ii)}=pout{ii};
            probs05{avls.indCurve(ii)}=p05out{ii};
            probs95{avls.indCurve(ii)}=p95out{ii};
        end
        load /cardinal3/r39g9/datsum/sumdataPOINTS.mat
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
                xvls=[times{ii}(crdayind) times{ii}(crdayind(end:-1:1))]-times{1}(1);
                yvls=[probs05{ii}(crdayind),probs95{ii}(crdayind(end:-1:1))];
                fill(xvls,yvls,'b')
                plot(times{ii}(crdayind)-times{1}(1),probs{ii}(crdayind),'r')
                end
            end   
        % Plot learning points + CIs
        for ii=avls.indPoint
                xvls=[times{ii};times{ii}(end:-1:1)]-times{1}(1);
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
        plot([t1-1-times{1}(1) t2+1-times{1}(1)],[p1 p1],'k')
        xlim([-0.5 4])
        ylim([0.3 0.9])
        
