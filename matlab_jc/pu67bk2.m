% pu67bk2


% syntax candidate
% Screened by TW
    % 05.18.10 --- surgery - bilateral RA cannulae implantation
% post-surgery song looks good

    % 05.20.10 --- probes in
% Branch = two downsweeps
% A (70%) - noisy note
% B (30%) - high stack


% 06.29
%       2mM apv on at 5:28pm at 1.5uL/min
%                  to 1.0uL/min at 5:45pm
%       acsf on at 7pm
% 06.30
%       2mM apv on at 5:50pm at 1.5uL/min
%      started going to shit around 7:15pm but got lots of data beforehand
%      accidentally turned WN off at 7:40pm
%       acsf on at 8:10pm
%      wn off for good at 8:10pm - intentional to allow recovery
% 07.01 - hanging upside down by toenails at 9:30am - released - JW to trim
        % trimmed at 10:30 am - gave him meloxicam too\
     % 0701pm and 0702am - definitely singing but rate seems to be lower
     % than before

% July 22, 2010 - analysis
% SYNSHIFT
% /syn_627_wnon
i=8;
%     Experiment1(i).listchars=notestats;
%         Experiment1(i).fvA=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
%     %Experiment1(i).fvB=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
%     Experiment1(i).fvD=findwnoteJC('batchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
%     Experiment1(i).fvE=findwnoteJC('batchnotes','e','','',0,[2000 2700],8500,1,'obs0',1);
%     Experiment1(i).fvC=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);

% Run fvals
%%%
    Experiment1(i).listcharscatch=notestats;
    Experiment1(i).fvAcatch=findwnoteJC('batchcatchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvBcatch=findwnoteJC('batchcatchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvDcatch=findwnoteJC('batchcatchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvEcatch=findwnoteJC('batchcatchnotes','e','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvCcatch=findwnoteJC('batchcatchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
% Learning pre
        pBaseA=0.952; % from experiment 1
    % apv pre
        pBaseAapv=0.921; % from experiment 2
  %%% Learning
    timesA=[];timesB=[];timesC=[];listchars=[];timechars=[];
    for i=[3 5 7] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesC=[timesC timing3(Experiment1(i).fvCcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='c'))];
    end
    LClearn=listchars;
    countA=0;countB=0;countC=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
                    end;end;end;
    end
    TClearn=timechars;
% Learning APV - second experiment - hit B
    for i=[4 6 8] % learning ACSF
        timesA=[];timesB=[];timesC=[];listchars=[];timechars=[];
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesC=[timesC timing3(Experiment1(i).fvCcatch)];
        LClearn_apv(i).data=Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='c'));
        listchars=LClearn_apv(i).data;
        countA=0;countB=0;countC=0;
        for j=1:length(listchars)
            if listchars(j)=='a';countA=countA+1;timechars(j)=timesA(countA);
            else if listchars(j)=='b';countB=countB+1;timechars(j)=timesB(countB);
                else if listchars(j)=='c';countC=countC+1;timechars(j)=timesC(countC);
                        end;end;end;
        end
        TClearn_apv(i).data=timechars;
    end

 %%%% 
 runanalysis(1-(LClearn=='a'),1,1-pBaseA,0.005,0)
 
         figure;hold on;
        % Experiment 1 - hit D, starting at baseline
            plot(TClearn,p05(1:end-1),'.')
            plot(TClearn,p95(1:end-1),'.')
            plot(TClearn,pmid(1:end-1))
            plot([min(TClearn)-20 max(TClearn)+20],[1-pBaseA 1-pBaseA],'k')
            for i=[4 6 8]
            plot(mean(TClearn_apv(i).data),mean(1-sum(LClearn_apv(i).data=='a')/length(LClearn_apv(i).data)),'+','Markersize',15,'Color','r')
            end
            plot(min(TClearn)-5,1-pBaseDapv,'+','Markersize',15,'Color','r')
            ylim([0 1])
            xlim([8255 8325])

