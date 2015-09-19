function pitch=findescapes(batch,batchnotes,low,high)
fvals=findwnoteJC(batchnotes,'a','','',0,[low high],8000,1,'obs0',1);
[vals,trigs]=triglabel(batch,'a',1,1,0,0);
count=0;
for i=1:length(trigs)
    for j=1:length(trigs(i).ntmintmp)
        count=count+1;
        if isempty(find(trigs(i).trigmintmp==trigs(i).ntmintmp(j)))
            miss(count)=1;
        else
            miss(count)=0;
        end
    end
end
ind=find(miss==1);
fvalsmiss=fvals(ind);
for i=1:length(fvalsmiss)
    shifted(i,:)=fvalsmiss(i).datt;
end
pitch=jc_pitchmat1024(shifted,1024,1020,2,low,high,1,'obs0',1);

            
