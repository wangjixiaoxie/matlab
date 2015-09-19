function [Datt]=getchandata(batchnotes,NOTES,CHANINDS,numdatapoints,pretimems)
% [Datt]=getchandata('batchnotes','ab',[0 2 4],1e4,50);
    %        [b1,a1]=butter(4,[30/32000],'high');
for i=1:length(NOTES)
    thisnote=NOTES(i);
    for j=1:length(CHANINDS)
        thischan=['obs' num2str(CHANINDS(j))];
        fvals=findwnoteSPK(batchnotes,thisnote,'','',0,[2000 2700],numdatapoints,1,thischan,1,pretimems);
        for k=1:length(fvals)
            Datt(i,j).data(k,:)=fvals(k).datt;%-mean(fvals(k).datt);
        end
    end
end
    