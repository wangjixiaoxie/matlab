pu58w25
% Reinforcement verification test - delayed WN - Feb. 2009
    % died after implant - March 2009
% Baseline data - February 1
% Experiment 1 - February 2 and first half of February 3
    % trig off the end of A, play ~30ms WN (goal was 13ms but smeared) with 22ms silence beforehand
    % A shifts, next note does not

    
load /cardinal3/pu58w25/Exp1data.mat
    figure;hold on
    clim=[7,16];
    colormap(hot)
    imagesc(t,f,log(avA),clim);syn;ylim([0,1e4]);xlim([-0.1 0.2])
    plot([1:1200]/8000-512/32000,mean(pitch203([1:1200],:)'),'Linewidth',4)
    %hold on;plot([950:1080],mean(pitch203([950:1080],:)')-mean(pitch201base([950:1080],:)'))
    hold on;plot([min(toffsall) max(toffsall)+240]/8000-512/32000,[3000 3000],'g','Linewidth',4)



    
% Experiment 2 - second half of February 3
    % trig off the end of A, play ~30ms WN with 32ms silence beforehand
    % 
% Experiment 3 - go back to 22ms silence