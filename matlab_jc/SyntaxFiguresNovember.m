clear all
load /bulbul3/AnalysisSyntax/birds.mat

%%% SYNTAX REVERSION SUMMARY DATA
                % for i=1:length(bird)
                %     baseline=bird(i).acsfmean(1);
                %     for j=2:length(bird(i).acsfmean)
                %         learnedAC(i,j-1)=bird(i).acsfmean(j)-baseline;
                %         learnedAPV(i,j-1)=bird(i).apvmean(j)-baseline;
                %         apvoff(i,j-1)=bird(i).apvmean(1)-bird(i).acsfmean(1);
                %         day(i,j-1)=bird(i).days(j);
                %     end
                % end
                % % apvoff (i.e. baseline apv offset) is n.s p=0.3;
                % learnedAPVadj=learnedAPV-apvoff;
                %         figure;hold on;
                %         for i=1:3
                %             plot(i-0.2,learnedAC(day==i),'k.','Markersize',15)
                %             plot(i+0.2,learnedAPVadj(day==i),'r.','Markersize',15)
                %             plot([(i-0.2) (i+0.2)],[learnedAC(day==i) learnedAPVadj(day==i)],'k')
                %         end
                %         figure;hold on;
                %         plot(learnedAC(learnedAC~=0),learnedAPVadj(learnedAC~=0),'r.','Markersize',15)
                %         plot([0 0.6],[0 0.6],'k')
                %         xlim([0 0.6])
                %         ylim([0 0.6])

                % Best figure
                    figure;hold on;
                    plot(learnedAC(day>1),learnedAPVadj(day>1),'r.','Markersize',15)
                    %plot(learnedAC(day==1),learnedAPVadj(day==1),'b.','Markersize',15)    
                    plot([0 0.7],[0 0.7],'k')
                    xlim([0 0.7])
                    ylim([0 0.7])
                  % n=15 inactivations in n=6 birds
                  % For days 2-3, we get n=11 inactivations in n=6 birds
% CV reduction for Syntax experiments
        clear all
        load /cardinal3/SyntaxLearningFigs/CV.mat
        
        
%%% FF REVERSION SUMMARY DATA    
clear all
load /cardinal3/SyntaxLearningFigs/birds.mat
load /cardinal3/SyntaxLearningFigs/birdsFF.mat % eventually put it all here
delAC=[];
delAPV=[];
count=1;
for i=1:length(birdFF)
    bF=birdFF(i);
    coef=1-2*(isequal(bF.direction,'down'));
    for j=2:length(bF.days)
        delAC(count)=coef*(bF.acsfmean(j)-bF.acsfmean(1));
        delAPV(count)=coef*(bF.apvmean(j)-bF.acsfmean(1));        
        count=count+1;
    end
end
%
delACSyn=[];
delAPVSyn=[];
count=1;
for i=1:length(bird)
    bS=bird(i);
    coef=1;
    for j=2:length(bS.days)
        delACSyn(count)=coef*(bS.acsfmean(j)-bS.acsfmean(1));
        delAPVSyn(count)=coef*(bS.apvmean(j)-bS.acsfmean(1));        
        count=count+1;
    end
end
figure;hold on;
subplot(121);hold on;bar([mean(delAC) mean(delAPV)]);
plot([1 1],[mean(delAC)-std(delAC)/length(delAC) mean(delAC)+std(delAC)/length(delAC)],'-')
plot([2 2],[mean(delAPV)-std(delAPV)/length(delAPV) mean(delAPV)+std(delAPV)/length(delAPV)],'-')
subplot(122);hold on;bar([mean(delACSyn) mean(delAPVSyn)])
plot([1 1],[mean(delACSyn)-std(delACSyn)/length(delACSyn) mean(delACSyn)+std(delACSyn)/length(delACSyn)],'-')
plot([2 2],[mean(delAPVSyn)-std(delAPVSyn)/length(delAPVSyn) mean(delAPVSyn)+std(delAPVSyn)/length(delAPVSyn)],'-')



%%%% Syntax example - r39g39 - see also SyntaxFiguresAPV2.m
        clear all
      % edit SyntaxAPVdata.m
        load /cardinal3/r39g9/datsum/avls.mat
        load /cardinal3/r39g9/datsum/sumdataCURVES.mat
        clear times probs probs05 probs95
        for ii=1:length(tms)
            times{avls.indCurve(ii)}=tms{ii};
            probs{avls.indCurve(ii)}=1-pout{ii};
            probs05{avls.indCurve(ii)}=1-p05out{ii};
            probs95{avls.indCurve(ii)}=1-p95out{ii};
        end
        load /cardinal3/r39g9/datsum/sumdataPOINTS.mat
        for ii=1:length(tms)
            times{avls.indPoint(ii)}=tms{ii};
            probs{avls.indPoint(ii)}=1-pout{ii};
            probs05{avls.indPoint(ii)}=1-p05out{ii};
            probs95{avls.indPoint(ii)}=1-p95out{ii};
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
        if length(times{1})==1;t1=times{1};else t1=min(times{1});end
        if length(times{end})==1;t2=times{end};else t2=max(times{end});end
        %for i=1:length(avls.indPointAC); pp(i)=probs{avls.indPointAC(i)}; end
        %if length(pp)==1;p1=pp;else p1=mean(pp);end
        plot([t1-1-times{1}(1) t2+1-times{1}(1)],[1-avls.pbase 1-avls.pbase],'k')
        xlim([-0.5 4])
        ylim([0 0.7])
        

%%%% FF example - r39g39 - see also SyntaxPitchFiguresAPV2.m
    clear all
            load /cardinal3/SyntaxBirds/FFsyntax1107.mat
            %%%% r39g39 - FF - hamming window
            datas=FFsyntax(3).data;
            for i=1:length(datas)
                acsfbase(i)=(datas(i).acsf & datas(i).baseline);
                acsfwn(i)=(datas(i).acsf & ~datas(i).baseline);
                apv(i)=(~datas(i).acsf);
            end
            acsfwnpitch=[];
            acsfwntimes=[];
            for i=1:length(datas)
                if acsfwn(i) | acsfbase(i);
                    times=timing3(datas(i).fvals);
                    pts=mean(datas(i).pitch(datas(1).window,:));
                    acsfwntimes=[acsfwntimes times];
                    acsfwnpitch=[acsfwnpitch pts];
                end
            end
            % Use Kalman filter to generate estimate of learning state and
            % confidence intervals.
                acjk=[zeros(1,50) acsfwnpitch-median(acsfwnpitch(1:50))]; % pre-processing
                [mn,se]=jckalman(acjk); % smooth
                mn=mn(51:end)+median(acsfwnpitch(1:50)); % post-processing
                se=se(51:end);
          
            figure;hold on;
            unidays=find(diff(acsfwntimes)>8);
            unistart=[1 unidays+1];
            uniend=[unidays length(acsfwntimes)];
            for i=1:length(unistart)
                indday=unistart(i):uniend(i);
                times=acsfwntimes(indday);
                mnday=mn(indday);
                selow=mn(indday)-2*se(indday);
                sehigh=mn(indday)+2*se(indday);
                unqtimes=unique(times);   
                unqselow=[];unqsehigh=[];unqmn=[];
                for j=1:length(unqtimes)
                    selind=find(times==unqtimes(j));
                    unqselow(j)=mean(selow(selind));
                    unqsehigh(j)=mean(sehigh(selind));
                    unqmn(j)=mean(mnday(selind));
                end
                fill([unqtimes unqtimes(end:-1:1)],[unqselow unqsehigh(end:-1:1)],'b')
                plot(unqtimes,unqmn,'r-')
            end
            for i=find(acsfbase+apv)
                times=mean(timing3(datas(i).fvals));
                pts=mean(datas(i).pitch(datas(1).window,:));
                ptmn(i)=mean(mean(datas(i).pitch(datas(1).window,:)));
                guess=[];
                for ii=1:1000
                    randsam=ceil(length(pts)*rand(1,length(pts)));
                    guess(ii)=mean(pts(randsam));
                end
                pt05=prctile(guess,5);
                pt95=prctile(guess,95);
                if ~acsfbase(i)
                    plot(times,ptmn(i),'r.','Markersize',15)
                    plot([times times],[pt05 pt95],'r-','Linewidth',2)
                end
            end
            plot([6080 6200],[mean([ptmn(1) ptmn(3)]) mean([ptmn(1) ptmn(3)])],'k')
            
                        
            
%%%% Syntax example 2 - r37g7
clear all
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
%%%% FF example 2 - r37g7
clear all
    load /cardinal3/SyntaxBirds/FFsyntax1107.mat
        %%%% r37g7 - FF - hamming window
            datas=FFsyntax(2).data;
            lenham=11;
            h=hamming(lenham);
            h=h*(lenham/sum(h));
            for i=1:length(datas)
                acsfbase(i)=(datas(i).acsf & datas(i).baseline);
                acsfwn(i)=(datas(i).acsf & ~datas(i).baseline);
                apv(i)=(~datas(i).acsf);
            end
            acsfwnpitch=[];
            acsfwntimes=[];
            for i=1:length(datas)
                if acsfwn(i);
                    times=timing3(datas(i).fvals);
                    pts=mean(datas(i).pitch(datas(1).window,:));
                    acsfwntimes=[acsfwntimes times];
                    acsfwnpitch=[acsfwnpitch pts];
                end
            end
            figure;hold on;
            unidays=find(diff(acsfwntimes)>8);
            unistart=[1 unidays+1];
            uniend=[unidays length(acsfwntimes)];
            for i=1:length(unistart)
                times=acsfwntimes(unistart(i):uniend(i));
                pts=acsfwnpitch(unistart(i):uniend(i));
                hamfilttimes=[];hamfiltpts=[];
                for ii=1:length(times)-lenham
                    hamfilttimes(ii)=mean(times(ii:ii+lenham-1).*h');
                    hamfiltpts(ii)=mean(pts(ii:ii+lenham-1).*h');
                end
                plot(hamfilttimes,hamfiltpts,'k.','Markersize',15)
            end
            for i=find(acsfbase+apv)
                times=mean(timing3(datas(i).fvals));
                pts=mean(datas(i).pitch(datas(1).window,:));
                ptmn(i)=mean(mean(datas(i).pitch(datas(1).window,:)));
                guess=[];
                for ii=1:1000
                    randsam=ceil(length(pts)*rand(1,length(pts)));
                    guess(ii)=mean(pts(randsam));
                end
                pt05=prctile(guess,5);
                pt95=prctile(guess,95);
            %     times=timing3(datas(i).fvals);
            %     pts=mean(datas(i).pitch(datas(1).window,:));
            %     hamfilttimes=[];hamfiltpts=[];
            %     for ii=1:length(datas(i).fvals)-lenham
            %         hamfilttimes(ii)=mean(times(ii:ii+lenham-1).*h');
            %         hamfiltpts(ii)=mean(pts(ii:ii+lenham-1).*h');
            %     end
                if acsfbase(i)
                    plot(times,ptmn(i),'b.','Markersize',15)
                    plot([times times],[pt05 pt95],'b-','Linewidth',2)
                else
                    plot(times,ptmn(i),'r.','Markersize',15)
                    plot([times times],[pt05 pt95],'r-','Linewidth',2)
                end
            end
                plot(hamfilttimes(55:end),hamfiltpts(55:end),'k.','Markersize',15)
            plot([2680 2830],[mean([ptmn(1) ptmn(3)]) mean([ptmn(1) ptmn(3)])],'k')
            xlim([2690 2830])