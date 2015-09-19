function [mmean,postesc,posthit]=jcTSA919(batch,pitch,fvals,back)
%Test this is the near future on the 917amTSANALY dataset
%[postesc,posthit]=jcTSA919('batchTSAnotcatch','batchTSAnotes',templa(:,1),cntrng(1));
%fvals=findwnote4(batchnotes,'a','','',0,[2000 2700],8500,1,'obs0',0);
fs=32000;
%mk_tempf(batch,templa,2,'obs0');
%get_trigt2(batch,cntrng,0.3,128,0,1);
%for i=1:length(fvals)
%    shifted(i,:)=fvals(i).datt;
%end
%[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,1950,2600,[1 ],'obs0',1);
[vals,trigs]=triglabel(batch,'a',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end





toffset=((toff/1000)*(fs)-512)/4;
count=0;

for i=1:length(fvals)
    if max(fvals(i).datt>5000)
        count=count+1;
        offset_time(i)=700; %round(toffset(count));
        hit(i)=1;
        
    else
        offset_time(i)=700; %round(toffset(round(1+rand*length(toffset))));
    end
    mmean(i)=median(pitch(offset_time(i)-200:offset_time(i)-1,i));
end
countesc=0;
counthit=0;
for j=1:length(hit)-back
    if hit(j)==0;
        countesc=countesc+1;
        postesc(countesc)=mmean(j+back);
    else
        counthit=counthit+1;
        posthit(counthit)=mmean(j+back);
    end
end
g=8;

postesc=median(postesc);
posthit=median(posthit);