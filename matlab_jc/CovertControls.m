    % /bulbul2/CovertAnalysis
    % /bulbul4/CovertAnalysis
    

%% CV % December 9, 2011

a=mean(std(Exps1(4).pitchACpre([180:300 900:1060],:)))/mean(mean(Exps1(4).pitchACpre([180:300 900:1060],:)))

b=mean(std(Exps1(4).pitchAPV([180:300 900:1060],:)))/mean(mean(Exps1(4).pitchAPV([180:300 900:1060],:)))
VarPre(4)=a;
VarAPV(4)=b;

varrecov=VarPost2(ind)./VarPre(ind)
    varreduc=VarAPV(ind)./VarPre(ind);
    figure;hold on;
    bar(0.5,1)
    bar(1.5,mean(varreduc))
    plot([1.5 1.5],[mean(varreduc)+std(varreduc)/sqrt(length(varreduc)) mean(varreduc)-std(varreduc)/sqrt(length(varreduc))],'Linewidth',2)
    plot(1.5,varreduc,'o','markersize',15)
    bar(2.5,mean(varrecov))
    plot([2.5 2.5],[mean(varrecov)+std(varrecov)/sqrt(length(varrecov)) mean(varrecov)-std(varrecov)/sqrt(length(varrecov))],'Linewidth',2)
    plot(2.5,varrecov,'o','markersize',15)
    ylim([0 1.5])
    

%%
    
% November 17, 2011
% AFP block - is AFP output completely blocked in experiments that showed
% learning?
load /bulbul4/CovertAnalysis/CovertCompact4.mat
    % Test 1: Does learning only occur in experiments with low variation
    % reduction?
            % Variation reduction
                varreduc=1-VarAPV(ind)./VarPre(ind);
            % Expression of Learning
                learnEXP=100*finalpoint(ind); % starts at 0.20%
            % Acquistion of Learning
                learnACQ=100*post2(ind); % starts at 1.0%

            [r,p]=corrcoef(learnACQ,varreduc) % R=-0.02
            figure;plot(varreduc,learnACQ,'.','Markersize',15)
            mean(learnACQ(find(varreduc>0.4))) % 1.18% - even better!!!
            [h,p]=ttest((learnACQ(find(varreduc>0.4))))   
                        % significant, p=0.002
            [h,p]=signtest((learnACQ(find(varreduc>0.4))))   
                        % significant, p=0.002
                        
            [h,p]=corrcoef(learnEXP,varreduc) % R=0.04
            figure;plot(varreduc,learnEXP,'.')
            mean(learnEXP(find(varreduc>0.4))) % 0.30% - partial expression not due to lack of var reduction
            [h,p]=ttest((learnEXP(find(varreduc>0.4))))
                        % not significant, p=0.20, n=10
            [h,p]=signtest((learnEXP(find(varreduc>0.4))))
                        % not significant, p=0.34, n=10
    % Test 2: Does learning only occur in experiments with expression?
    
        [h,p]=corrcoef(learnEXP,learnACQ) % R=0.19
        figure;plot(learnEXP,learnACQ,'.','Markersize',15)
        hold on;plot([-1.5 1.5],[0 0])
        plot([0 0],[-0.5 3.5])
        xlim([-1.5 1.5]);ylim([-0.5 3.5])
        
        % note that the following is a very harsh control because the average expression
        % in these cases is -0.49% - not at all centered around zero
            % and yet it still works out!!!
  
        mean(learnACQ(find(learnEXP<0))) % 0.87% - pretty good, n=8
        figure;hold on;
        bar(1,mean(learnACQ))
        bar(2,mean(learnEXP))
        plot([1 1],[mean(learnACQ)-std(learnACQ)/sqrt(21) mean(learnACQ)+std(learnACQ)/sqrt(21)],'-','Linewidth',2)        
        plot([2 2],[mean(learnEXP)-std(learnEXP)/sqrt(21) mean(learnEXP)+std(learnEXP)/sqrt(21)],'-','Linewidth',2)
        bar(4,mean(learnACQ(find(learnEXP<0))))
        plot([4 4],[mean(learnACQ(find(learnEXP<0)))-std(learnACQ(find(learnEXP<0)))/sqrt(8) mean(learnACQ(find(learnEXP<0)))+std(learnACQ(find(learnEXP<0)))/sqrt(8)],'-','Linewidth',2)
        bar(5,mean(learnEXP(find(learnEXP<0))))
        plot([5 5],[mean(learnEXP(find(learnEXP<0)))-std(learnEXP(find(learnEXP<0)))/sqrt(8) mean(learnEXP(find(learnEXP<0)))+std(learnEXP(find(learnEXP<0)))/sqrt(8)],'-','Linewidth',2)
        xlim([0 5])
        %
        figure;hold on;
        bar(2,mean(learnACQ(find(learnEXP<0))))
        plot([2 2],[mean(learnACQ(find(learnEXP<0)))-std(learnACQ(find(learnEXP<0)))/sqrt(8) mean(learnACQ(find(learnEXP<0)))+std(learnACQ(find(learnEXP<0)))/sqrt(8)],'-','Linewidth',2)
        plot(2,(learnACQ(find(learnEXP<0))),'.','markersize',15)
        bar(0.5,mean(learnEXP(find(learnEXP<0))))
        plot([0.5 0.5],[mean(learnEXP(find(learnEXP<0)))-std(learnEXP(find(learnEXP<0)))/sqrt(8) mean(learnEXP(find(learnEXP<0)))+std(learnEXP(find(learnEXP<0)))/sqrt(8)],'-','Linewidth',2)
        plot(0.5,(learnEXP(find(learnEXP<0))),'.','markersize',15)
        xlim([0 3])

        [h,p]=ttest((learnACQ(find(learnEXP>0))))
                        %  significant, p<0.001, n=8
        [h,p]=signtest((learnACQ(find(learnEXP>0))))
                        %  significant, p=0.0034, n=8
        
    % Test 3: Can any linear combination of the following predict acquisition?
        % - lack of variability reduction (1-varreduc)
        % - expression of learning (learnEXP)
        [b,bint,r,rint,stats] = [b,bint,r,rint,stats] = regress(learnACQ',[ones(1,21);1-varreduc;learnEXP]')
            % stats - note that for this to make since I need the column of ones as a regressor (as I did include above)
                % R^2-stat= 0.0363
                % F-stat = 0.3392
          ##    % p-value for F-stat=0.7168
                % estimate of error variance=0.6819
                
 % Test 4: look at cases with no variation and no expression? Too small
 % sample size?
                    