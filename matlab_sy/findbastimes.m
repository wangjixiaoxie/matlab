function [avls]=findbastimes(avls)
    
    dirvls=getvals(avls.fvcombdir{1},1,'trig');
    
    flvls=floor(dirvls(:,1));
    uniquetims=unique(flvls);
    avls.dirdays(:,1)=(uniquetims+7/24);
    avls.dirdays(:,2)=(uniquetims+21/24);