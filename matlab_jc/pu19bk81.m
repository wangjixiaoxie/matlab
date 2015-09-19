% pu19bk81
% MMAN lesion 08.26.11 (friday)


cleandir4('batch',1e4,500,5,5);
% started singing to female 08.30.11
% started singing alone 08.31.11
% Increased stereotypy at the branch point!!!

% 9.1.11 - 9.2.11 - 96% a-->b, 4% a--c or a--x
    % template needs work re: false negatives
% Looks like sequence learning but incomplete recovery (unstable baseline?)
    
% 9.14 --> adjusted template
% 9.15 at dawn - wn on - hit dominant transition
% 9.18 - wn off at dark

% 9.21 - noon - wn on - hit below 3230Hz
% 9.22 - 1pm - wn off

load DataFF0922.mat
% average off morning before wn (9.21) to morning at conclusion of wn (9.22)
mean(median(pitchBwn(180:280,124:end)')-median(pitchBpre(180:280,371:end)'))
figure;plot(timing3(fvalsBpre),mean(pitchBpre(200:250,:)),'.')
hold on;plot(timing3(fvalsBwn),mean(pitchBwn(200:250,:)),'r.')
hold on;plot(timing3(fvPost),mean(pitchPost(200:250,:)),'k.')


load Data0919.mat
    figure;hold on;
     plot(runningaverage(tvsorted912,50),runningaverage(numsorted912,50),'Color','b')
     plot(runningaverage(tvsorted914,50),runningaverage(numsorted914,50),'Color','b')
   
    plot(runningaverage(tvsorted915,50),runningaverage(numsorted915,50),'Color','r')
    plot(runningaverage(tvsorted920,50),runningaverage(numsorted920,50),'Color','k')

dirf('*.cbin.not.mat','batchnotes')
fvalsB=findwnoteJC('batchnotes','b','a','',0,[2000 2700],8500,1,'obs0',1)
fvalsC=findwnoteJC('batchnotes','c','a','',0,[2000 2700],8500,1,'obs0',1)
%fvalsX=findwnoteJC('batchnotes','x','a','',0,[2000 2700],8500,1,'obs0',1)
tv914B=timing4(fvalsB);
tv914C=timing4(fvalsC);
%tv914X=timing4(fvalsX);
tvall914=[tv914B tv914C ];
numsall2=[ones(1,length([tv914B])) zeros(1,length([tv914C])) ];
[b2,ind2]=sort(tvall914);
numsorted914=numsall2(ind2);
tvsorted914=tvall914(ind2);

