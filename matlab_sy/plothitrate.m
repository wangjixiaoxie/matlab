%plothitrate.m

%for each day plot the mean hit rate
%and plot the standard deviation.

initday=datenum(2006,11,21,11,54,0)

init=datenum(2006,11,16,11,0,0)


vals=getvals(fv,2,'TRIG');
ind=find(vals(:,1)>initday)

valsrstr=fix(vals(ind,:));

dayind=unique(valsrstr(:,1))

for i=1:length(dayind)
   indtrig=find(valsrstr(:,1)==dayind(i)&valsrstr(:,3)==1)
   indnotrig=find(valsrstr(:,1)==dayind(i)&valsrstr(:,3)==0)
   hitrate(i)=length(indtrig)/(length(indtrig)+length(indnotrig))
end
figure

plot(dayind-fix(init)+1,hitrate,'.')


