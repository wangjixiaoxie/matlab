function vals=get_songlengths(batch,NT);%,MAXINT);
%vals=get_songlengths(batch,NT,MAXINT);
%batch=batchfile
%NT=note to use
%MAXINT in seconds intervals bigger than this = newsong


fid=fopen(batch,'r');
fns=[];cnt=0;
while (1)

    fn=fgetl(fid);

    if (~ischar(fn))
        break;
    end

    if (~exist([fn,'.not.mat'],'file'))
        continue;
    end
    cnt=cnt+1;
    fns(cnt).name=fn;
end
fclose(fid);
vals=[];
for ii=1:length(fns)
    fn=fns(ii).name;
    load([fn,'.not.mat']);

    pp=findstr(labels,NT);
    if (length(pp)>1)
        dt=diff(onsets(pp))*1e-3;

        vals=[vals;length(pp),(offsets(pp(end))-onsets(pp(1)))*1e-3];
    end
end
return;