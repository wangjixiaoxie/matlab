function []=plotlabeledspect(fn,timebnds,freqbnds,label)
    [dat,fs]=ReadCbinFile(fn);
    evspect(dat(:,1),fs,[ freqbnds(1) freqbnds(2)]);
    hold on;
    notmatfn=[fn '.not.mat'];
    cmd=['load ' notmatfn];
    eval(cmd);
    
    if(label)
    %onsets and offsets are in milliseconds
    for ii=1:length(onsets)
        text(((onsets(ii)+offsets(ii))/2000),1000,labels(ii),'Color','r','FontSize',18);
    end
    end