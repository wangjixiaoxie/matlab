
% FOUND DEAD 3.28.11


w35pk25
% NO SYNTAX (no branch points)
% DECENT FOR COVERT - one note, maybe two - positive control done

% refrac= 0.2

% segment at 500

% B low stacks followed by A high stacks which are sometimes noisy

%%%%% Experiment 1
% Hit B low stacks EARLY @ 100% 
% 04.07.10 - template made - BTAF off of downsweep before B
% 
% 04.09.10 - hit 100% with catch trials (above 2250Hz)
    % wn on at 8pm on 04.08
% 04.12.10 - power was out twice in the morning, for a total of 3hrs, but
% hits are still good
    % labeled 50% of catch trials'
    % 2550 - 9pm
    % 2560 - 11am
 % Allow learning for several days - needs to be a big shift to use as the 
 % example figure 
figure;plot(tv407B,pitch407B(200,:),'*')
hold on;plot(tv408B,pitch408B(200,:),'*','Color','k')
plot(tv411wn,pitch411wn(200,:),'*','Color','r')
plot(tv413wn,pitch413wn(200,:),'*','Color','r')
plot(tv415,pitch415(200,:),'*')
plot(tv418A,pitch418A(200,:),'*')

window=[160:220];
figure;plot(runningaverage(tv407B,50),runningaverage(mean(pitch407B(window,:)),50),'*')
hold on;plot(runningaverage(tv408B,50),runningaverage(mean(pitch408B(window,:)),50),'*','Color','k')
plot(runningaverage(tv411wn,50),runningaverage(mean(pitch411wn(window,:)),50),'*','Color','r')
plot(runningaverage(tv413wn,50),runningaverage(mean(pitch413wn(window,:)),50),'*','Color','r')
plot(runningaverage(tv415,50),runningaverage(mean(pitch415(window,:)),50),'*')
plot(runningaverage(tv418A,50),runningaverage(mean(pitch418A(window,:)),50),'*')

% 07.09.10 - implanted bilateral cannulae in LMAN
% 07.12.10 - inserted probes at 12:30pm with Ket/Mid anesthetic
% 07.13.10 - 200um muscimol on at 10:25am at 1.5uL/min
    % acsf on at 12:20pm at 1.5uL/min
    % 45% reduction
    mean(std(pitch713B(200:400,end-50:end)'))/mean(mean(pitch713B(200:400,end-50:end)'))
    mean(std(pitch713apvB(200:400,40:end)'))/mean(mean(pitch713apvB(200:400,40:end)'))
    mean(std(pitch713B2(200:400,:)'))/mean(mean(pitch713B2(200:400,:)'))
 figure;plot(timing3(fv713B),mean(pitch713B(200:400,:)),'*')   
 hold on;plot(timing3(fv713apvB),mean(pitch713apvB(200:400,:)),'*','Color','r')
 plot(timing3(fv713B2),mean(pitch713B2(200:400,:)),'*')   
% probes out
% 07.19.10 - probes in at 6pm
% 07.20.10 - 12:45pm - 200um muscimol on

% 07.27.10 - probes in at 7pm - independent lines
% 07.28.10 - cleaned probe on right side and re-inserted at 10am