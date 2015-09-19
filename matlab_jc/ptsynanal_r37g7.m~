
                        % clear fv fbins bt smoothamp baspath conbins avls
                        % avls.baspath='/cardinal4/SyntaxAPV/r37g7/'
                        % 
                        % % Syntax data 
                        % load /cardinal4/SyntaxAPV/r37g7/ExperimentFiles.mat
                        % inds1=[5 7 4 8 9 10];
                        % for i=1:6
                        %     avls.pvls{i}=Experiment1(inds1(i)).folder;  %%%
                        %     avls.cvl{i}='batchfiles'; 
                        % end
                        %     avls.WINDOWSIZE=16;
                        %     
                        %     
                        % % Which notes are hit with WN to reduce their probability?
                        %     avls.targnts(1)='a';
                        %     avls.NT{1}='a';
                        %     avls.NT{2}='b';
                        %     avls.NT{3}='d';
                        %     avls.NT{4}='e'
                        % 
                        % %%%%%%%
                        % avls.indCurve=[4 6];
                        % avls.indCurveAC=1:length(avls.indCurve);
                        % avls.indCurveAPV=[];
                        % avls.indPoint=[1:3 5];
                        % avls.indPointAC=[1 3];
                        % avls.indPointAPV=[2 4];
                        % 
                        % 
                        % % 3. Analyze and plot data
                        %         % Calculate learning curves
                                    avls1=avls;
                                    avls1.pvls=avls.pvls(avls.indCurve);
                                    avls1.cvl=avls.cvl(avls.indCurve);
                                    avls1.MKSYN=1:length(avls1.pvls);
                                    syn_settimestampslabelsCURVE(avls1) % saved as /datsum/sumdataCURVES.mat
                        %         % Calculate learning points + CIs
                        %             avls2=avls;
                        %             avls2.pvls=avls.pvls(avls.indPoint);
                        %             avls2.cvl=avls.cvl(avls.indPoint);
                        %             avls2.MKSYN=1:length(avls2.pvls);    
                        %             syn_settimestampslabelsPOINT(avls2) % saved as /datsum/sumdataPOINTS.mat
load /cardinal4/SyntaxAPV/r37g7/datsum/avls.mat   
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
                xvls=[times{ii}(crdayind);times{ii}(crdayind(end:-1:1))]-times{1};
                yvls=[probs05{ii}(crdayind) probs95{ii}(crdayind(end:-1:1))];
                fill(xvls,yvls,'b')
                plot(times{ii}(crdayind)-times{1},probs{ii}(crdayind),'g','Linewidth',2)
                end
            end   
        % Plot learning points + CIs
        for ii=avls.indPoint
                xvls=[times{ii};times{ii}(end:-1:1)]-times{1};
                if ii==1;xvls=[22.5;22.5];end
                if ii==2;xvls=[23;23];end
                if ii==3;xvls=[23.5;23.5];end
                yvls=[probs05{ii} probs95{ii}(end:-1:1)];
                if isempty(find(avls.indPointAC==ii))
                    plot(xvls,yvls,'r-','Linewidth',2)
                    plot(mean(xvls),mean(yvls),'r.','Markersize',15)
                else
                    plot(xvls,yvls,'b-','Linewidth',2)
                    plot(mean(xvls),mean(yvls),'b.','Markersize',15)
                end
        end

        for i=1:length(avls.indPointAC); pp(i)=probs{avls.indPointAC(i)}; end
        if length(pp)==1;p1=pp;else p1=mean(pp);end
        plot([22 26],[p1 p1],'k')
        xlim([22 26])