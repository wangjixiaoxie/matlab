tvallapvPRE=[tvalsapvPREA tvalsapvPREB];
numsall2=[zeros(1,length(tvalsapvPREA)) ones(1,length(tvalsapvPREB))];
[b2,ind2]=sort(tvallapvPRE);
numsortedapvPRE=numsall2(ind2);
tvsortedapvPRE=tvallapvPRE(ind2);
hold on;plot(runningaverage(tvallapvPRE,80),runningaverage(numsortedapvPRE,80),'*','Color','r')

tvallacsfPRE=[tvalsacsfPREA tvalsacsfPREB];
numsall2=[zeros(1,length(tvalsacsfPREA)) ones(1,length(tvalsacsfPREB))];
[b2,ind2]=sort(tvallacsfPRE);
numsortedacsfPRE=numsall2(ind2);
tvsortedacsfPRE=tvallacsfPRE(ind2);
figure;plot(runningaverage(tvallacsfPRE,80),runningaverage(numsortedacsfPRE,80),'*')

tvallacsf1=[tvalsacsf1A tvalsacsf1B];
numsall2=[zeros(1,length(tvalsacsf1A)) ones(1,length(tvalsacsf1B))];
[b2,ind2]=sort(tvallacsf1);
numsortedacsf1=numsall2(ind2);
tvsortedacsf1=tvallacsf1(ind2);
figure;plot(runningaverage(b2,80),runningaverage(numsortedacsf1,80),'*')

tvallacsf=[tvalsacsfA tvalsacsfB];
numsall2=[zeros(1,length(tvalsacsfA)) ones(1,length(tvalsacsfB))];
[b2,ind2]=sort(tvallacsf);
numsortedacsf=numsall2(ind2);
tvsortedacsf=tvallacsf(ind2);
hold on;plot(runningaverage(tvallacsf,80),runningaverage(numsortedacsf,80),'*')

tvallapv=[tvalsapvA tvalsapvB];
numsall2=[zeros(1,length(tvalsapvA)) ones(1,length(tvalsapvB))];
[b2,ind2]=sort(tvallapv);
numsortedapv=numsall2(ind2);
tvsortedapv=tvallapv(ind2);
hold on;plot(runningaverage(tvallapv,80),runningaverage(numsortedapv,80),'*','Color','r')


tvallrev=[tvalsrevA tvalsrevB];
numsall2=[zeros(1,length(tvalsrevA)) ones(1,length(tvalsrevB))];
[b2,ind2]=sort(tvallrev);
numsortedrev=numsall2(ind2);
tvsortedrev=tvallrev(ind2);
figure;plot(runningaverage(tvsortedrev,60),runningaverage(numsortedrev,60),'*','Color','r')

tvallapv2=[tvalsapv2A tvalsapv2B];
numsall2=[zeros(1,length(tvalsapv2A)) ones(1,length(tvalsapv2B))];
[b2,ind2]=sort(tvallapv2);
numsortedapv2=numsall2(ind2);
tvsortedapv2=tvallapv2(ind2);
figure;plot(runningaverage(tvsortedapv2,50),runningaverage(numsortedapv2,50),'*','Color','r')

tvallreverse2=[tvreverse2A tvreverse2B];
numsall2=[zeros(1,length(tvreverse2A)) ones(1,length(tvreverse2B))];
[b2,ind2]=sort(tvallreverse2);
numsortedreverse2=numsall2(ind2);
tvsortedreverse2=tvallreverse2(ind2);
figure;plot(runningaverage(tvsortedreverse2,50),runningaverage(numsortedreverse2,50),'*','Color','r')

tvallPRE1=[tvPRE1a tvPRE1c];
numsall2=[zeros(1,length(tvPRE1a)) ones(1,length(tvPRE1c))];
[b2,ind2]=sort(tvallPRE1);
numsortedPRE1=numsall2(ind2);
tvsortedPRE1=tvallPRE1(ind2);
figure;plot(runningaverage(tvsortedPRE1,50),runningaverage(numsortedPRE1,50),'*','Color','r')

tvall1103pre1=[tv1103pre1A tv1103pre1B];
numsall2=[zeros(1,length(tv1103pre1A)) ones(1,length(tv1103pre1B))];
[b2,ind2]=sort(tvall1103pre1);
numsorted1103pre1=numsall2(ind2);
tvsorted1103pre1=tvall1103pre1(ind2);
figure;plot(runningaverage(tvsortedPRE1,50),runningaverage(numsortedPRE1,50),'*','Color','r')

% This plot demonstrates the reduction in variability
        figure;subplot(211);hold on;
        plot(tvalsacsf1A,pitchacsf1A(250,:),'*')
        plot(tvalsapvA,pitchapvA(250,:),'*','Color','r')
        plot(tvalsacsfA,pitchacsfA(250,:),'*')




% SUMMARY PLOT
        figure;hold on;
        % Earliest acsf files (1103tmptestPM - previous manipulation of any form was >24hrs beforehand (apv on 11/01)
           plot(runningaverage(tvsorted1103pre1+30,30),runningaverage(numsorted1103pre1,30),'*','Color','b')
        % Initial apv run at baseline - to check for acute effect (1103_4mMapv)
            plot(runningaverage(tvsortedapvPRE+30,30),runningaverage(numsortedapvPRE,30),'*','Color','r')
        % Baseline - morning - just before WN turned on (1104tmptestJC)
            plot(runningaverage(tvsortedacsfPRE,30),runningaverage(numsortedacsfPRE,30),'*','Color','k')
        % White noise on with ACSF in (1105wnon)
            plot(runningaverage(tvsortedacsf1,30),runningaverage(numsortedacsf1,30),'*')
        % Inactivate with WN still on (1105_APVon)
            plot(runningaverage(tvsortedapv,30),runningaverage(numsortedapv,30),'*','Color','r')
        % ACSF back in with WN on (1105_acsfon)
            plot(runningaverage(tvsortedacsf(1:72),30),runningaverage(numsortedacsf(1:72),30),'*')
            plot(runningaverage(tvsortedacsf(73:end),30),runningaverage(numsortedacsf(73:end),30),'*')
        % Reverse the contingency (1106wn_revconting)
            plot(runningaverage(tvsortedreverse2(1:210),30),runningaverage(numsortedreverse2(1:210),30),'*','Color','k')
            plot(runningaverage(tvsortedreverse2(211:end),30),runningaverage(numsortedreverse2(211:end),30),'*','Color','k')
        % Inactivate with WN still on (1107_APVon)
            plot(runningaverage(tvsortedapv2,30),runningaverage(numsortedapv2,30),'*','Color','r')
        