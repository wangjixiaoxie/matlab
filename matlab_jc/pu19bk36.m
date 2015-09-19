% pu19bk36

% syntax/lesion candidate

% started screening 10.04.11
% decent song
% A-->B (high stacks) 36% of time
% A-->C (kind of noisy note) 64% of time
% template to hit C, has some false positives on a more noisy note
    % adjusted template at noon on 10.05 - looking good
    % wn on at dark on 10.05 (planned)
    
% 10.06 - wn on
    % 11:05am - adjusted template to avoid B
    % it's going in the wrong direction
    
    
 dirf('*.cbin.not.mat','batchnotes')
fvWNx=findwnoteJC('batchnotes','x','','',0,[2000 2700],8500,1,'obs0',1);
fvWNb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvWNc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvWNx=timing3(fvWNx);
tvWNb=timing3(fvWNb);
tvWNc=timing3(fvWNc);
tvallWN=timing3([fvWNx fvWNb fvWNc]);
numsall2=[zeros(1,length(tvWNx)) zeros(1,length(tvWNb)) ones(1,length(tvWNc))];
[b2,ind2]=sort(tvallWN);
numsortedWN=numsall2(ind2);
tvsortedWN=tvallWN(ind2);

fvPREx=findwnoteJC('batchnotes','x','','',0,[2000 2700],8500,1,'obs0',1);
fvPREb=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
fvPREc=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
tvPREx=timing3(fvPREx);
tvPREb=timing3(fvPREb);
tvPREc=timing3(fvPREc);
tvallPRE=timing3([fvPREx fvPREb fvPREc]);
numsall2=[zeros(1,length(tvPREx)) zeros(1,length(tvPREb)) ones(1,length(tvPREc))];
[b2,ind2]=sort(tvallPRE);
numsortedPRE=numsall2(ind2);
tvsortedPRE=tvallPRE(ind2);

figure;hold on;
plot(runningaverage(tvsortedPRE,60),runningaverage(numsortedPRE,60),'Color','b')
plot(runningaverage(tvsortedWN,60),runningaverage(numsortedWN,60),'Color','r')


   