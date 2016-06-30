% === compare the shift from single targ phase to the first two targ phase,
% for both congr and incongr experiments. Get the differences. compare to
% see if the difference for cong is more positive than difference for
% incongr

Xcongr=[-0.01, 0.5]; % bk34 and rd23 increases in consolidation
Yincongr=[-1. -0.15]; % same birds, but for incongr

[~,p]=ttest(Xcongr, Yincongr)
    
% === unapired, add all the other incongr experiments
Yincongr_unpaired=[Yincongr, -0.36, -0.07, -0.01, -0.2, 0.48];
[~, p]=ttest2(Xcongr, Yincongr_unpaired)

