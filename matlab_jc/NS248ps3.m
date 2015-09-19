%
% Problem 1
% A
% Simulate draws from the normal distns, as in last weeks PSET #5
sigma1=10;
sigma2=10;
mn1=20;
mn2=30;
r1=normrnd(mn1,sigma1,100000,1);
r2=normrnd(mn2,sigma2,100000,1);
for i=1:1000
    samps1(i,:)=r1(ceil(rand(1,10)*100000));
    samps2(i,:)=r2(ceil(rand(1,10)*100000));
end
mnsamps1=mean(samps1');
mnsamps2=mean(samps2');
sum(mnsamps1<mnsamps2)/length(mnsamps1);
%       98.9 percent of the time 

% B
clear all
sigma1=10;
sigma2=10;
mn1=25;
mn2=25;
r1=normrnd(mn1,sigma1,100000,1);
r2=normrnd(mn2,sigma2,100000,1);
for i=1:1000
    samps1(i,:)=r1(ceil(rand(1,10)*100000));
    samps2(i,:)=r2(ceil(rand(1,10)*100000));
end
for i=1:1000
    significant(i)=ttest2(samps1(i,:),samps2(i,:));
end
sum(significant)/length(significant);
%  alpha=0.05 ---     4.9 percent of the time
for i=1:1000
    significant(i)=ttest2(samps1(i,:),samps2(i,:),0.1);
end
%  alpha=0.1  ---     9.3 percent of the time
%%%% alpha is the false positive rate %%%%

% C
for i=1:1000
    significant(i)=ttest2(samps1(i,:),samps2(i,:),0.05);
end
% Choose four results, what is the probability that one or more will be
% one.
for i=1:250
    sumsign(i)=sum(significant(i*4-3:i*4))>=1;
end
%   sum(sumsign)/length(sumsign)=0.18           0.049*4=0.196
% run it again
%   sum(sumsign)/length(sumsign)=0.164          0.047*4=0.188

% Conclusion --- false positive probability is much larger than 0.05, but
% still a bit less than (false positive probability for a single test)*4

% Problem 2
% illustration of anovan function


% three groups, three separate times per group
g1{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g1{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g1{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

g2{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g2{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g2{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

g3{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g3{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
g3{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

y = [];
group = [];
time = [];
for i = 1:3
    y = [y g1{i}];
    group = [group ones(size(g1{i}))];
    time = [time i*ones(size(g1{i}))];
end
for i = 1:3
    y = [y g2{i}];
    group = [group 2*ones(size(g2{i}))];
    time = [time i*ones(size(g2{i}))];
end
for i = 1:3
    y = [y g3{i}];
    group = [group 3*ones(size(g3{i}))];
    time = [time i*ones(size(g3{i}))];
end

%%%% A.
yA=y;
yA(1:30)=yA(1:30)+10;
[p atab stats] = anovan(yA, {group time}, 3, 3, strvcat('group', 'time'));
comp = multcompare(stats, .05, 'off',[] ,[] ,[1 2])
% Analysis of Variance
%   Source	Sum Sq.	   d.f.	Mean Sq.	F	    Prob>F
%   group	   1781.3122	2	890.6561  93.6297	0
%   time	    9.7628	    2	4.8814	  0.51315	0.60054
%   group*time 65.4606	4	16.3652	  1.7204	0.15349
%   Error	    770.5155	81	9.5125		
%   Total	    2627.0511	89			


%%%% B.
yB=y;
yB(1:10)=yB(1:10)+10;
yB(31:40)=yB(31:40)+10;
yB(61:70)=yB(61:70)+10;
[p atab stats] = anovan(yB, {group time}, 3, 3, strvcat('group', 'time'));
comp = multcompare(stats, .05, 'off',[] ,[] ,[1 2])
%%%
% Analysis of Variance
%   Source	Sum Sq.	   d.f.	Mean Sq.	F	     Prob>F
%   group	    31.1595	    2	15.5798	  1.6378	0.20078
%   time	    1736.262	2	868.131	  91.2618	0
%   group*time 65.4606	4	16.3652	  1.7204	0.15349
%   Error	     770.5155	81	9.5125		
%   Total	    2603.3976	89			


%%%% C.
yC=y;
yC(31:40)=yC(31:40)+3;
[p atab stats] = anovan(yC, {group time}, 3, 3, strvcat('group', 'time'));
comp = multcompare(stats, .05, 'off',[] ,[] ,[1 2])
%%%
% Analysis of Variance
%   Source	Sum Sq.	   d.f.	Mean Sq.	F	    Prob>F
%   group	    26.2167	    2	13.1083	   1.378	0.25793
%   time	    2.4127	    2	1.2064	   0.12682	0.88107
%   group*time 113.988	4	28.497	   2.9957	0.023284
%   Error	    770.5155    81	9.5125		
%   Total	     913.133	89			


% Problem 3
% A. 
% Group 1 is increasing steadily (r=0.5, p=0.02), which makes things
% tricky.  Gaussian assumption is violated (lillietest).  
% Don't know if data are paired or not.
% Wilcoxon rank-sum test.
% B.
[p,h,k]=ranksum(group1,group2,0.05)
p=0.004  % YES
% C.
% Group 1 may have a memory deficit.  
% Regarding the increase: This deficit may get worse with age or may have a
% circadian relationship (assuming that values later in the group are tested 
% later in the day).  But note also that when you include data from both 
% groups, the trend is not significant (tested with ANOVA).

% Problem 4
% Group 1 is increasing again.  % Data do not look very
% Gaussian.
% Wilcoxon rank-sum test.
% B.
[p,h,k]=ranksum(group1,group2,0.05)
p=0.0155 % YES
% C.
% Same as Problem 3
% 

% Problem 5
% A. Same as above
% B. 
p=0.002 % YES
% C. Same as above

% Problem 6
% A.
totspikes1=zeros(1,length(spikes1));
totspikes2=zeros(1,length(spikes2));
for i=1:length(spikes1)
    % find all spikes in the first 500ms
    totspikes1(i)=sum(spikes1{i}<0.5);
    totspikes2(i)=sum(spikes2{i}<0.5);
end
% multiply by two to convert spikes in 500ms to Hz
% MEANS
rate1=mean(totspikes1*2); % 4.92Hz
rate2=mean(totspikes2*2); % 0.28Hz
% VARIANCES
var1=std(totspikes1*2); % 3.08Hz
var2=std(totspikes2*2); % 0.70Hz

% B.
sdest1=3.0828;
sdest2=0.7010;
mnest1=4.92;
mnest2=0.28;
Allspikes1=[];
Allspikes2=[];
for i=1:50
    Allspikes1=[Allspikes1 spikes1{i}];
    Allspikes2=[Allspikes2 spikes2{i}];
end
Allspikes1sort=sort(Allspikes1);
Allspikes2sort=sort(Allspikes2);
tot1(1)=find(Allspikes1sort<0.1,1,'last');
between1(1)=tot1(1);
tot2(1)=find(Allspikes2sort<0.1,1,'last');
between2(1)=tot2(1);
k=[0.1:0.1:2];
for i=2:length(k)
    times(i)=k(i);
    tot1(i)=find(Allspikes1sort<times(i),1,'last');
    between1(i)=tot1(i)-tot1(i-1);
    tot2(i)=find(Allspikes2sort<times(i),1,'last');
    between2(i)=tot2(i)-tot2(i-1);
end
zscore1=(0.02*between1*10-mnest1)/sdest1;
zscore2=(0.02*between2*10-mnest2)/sdest2;
% zscore1 =
%     0.1557   -0.1038   -0.1038   -0.2336    0.2855    
%     4.5024    3.3346    1.4532    1.3883    0.5450   
%     0.2206    0.1557   -0.0389    0.0260   -0.1038    
%     0.3503    0.4152    0.7396   -0.1038    0.3503

% zscore2 =
%     0.1712   -0.1141   -0.1141    0.1712   -0.1141   
%     32.1255   20.9986   13.5806   11.8688    9.3010    
%     7.0185    7.0185    7.5892    9.5863    8.7304    
%     6.4479    5.8773    5.8773    6.4479    8.7304

% C
% Raw data:
    % Cell 1 has 573 spikes between 0.5s and 2s
    % Cell 1 has 123 spikes during baseline

    % Cell 2 has 587 spikes between 0.5s and 2s
    % Cell 2 has only 7 spikes during baseline

% Both cells have similar activity after the stimulus, but Cell 2 has very
% low baseline activity, so it is more strongly modulated.

% Issues: use of a z-score assumes you know the population standard
% deviation, but you actually only know the sample standard deviation.
% Furthermore, if the data are not distributed in a Gaussian (baseline data
% from sample 2 clearly are not --- only 0Hz or 2Hz), then sample standard
% deviation will be a poor estimate.
