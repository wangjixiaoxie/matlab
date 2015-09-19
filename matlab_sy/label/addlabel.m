%outdated...current code is changelabel.
function trigs=addlabel(batch,exlabel,newlabel,post);


ff=load_batchf(batch);
for ifn=1:length(ff)
    fn=ff(ifn).name;
    fnn=[fn,'.not.mat'];
    if (~exist(fnn,'file'))
        continue;
    end
    disp(fn);
    load(fnn)
    labels = lower(labels);
    labels(findstr(labels,'0'))='-';
    [labelind]=findstr(labels,exlabel);
     if(post)
         labels(labelind+1)=newlabel;
     else
         labels(labelind-1)=newlabel;
     end
     cmd=['save -append ' fnn ' labels'];
     eval(cmd)
end