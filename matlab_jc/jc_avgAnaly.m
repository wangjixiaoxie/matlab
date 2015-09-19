function chisqLM=jc_avgAnaly(avLMAN,avNone,starting,ending,n)
chisqLM=0;
chisqNO=0;
for i=starting:ending
    chisqLM=chisqLM+abs((avLMAN(i)-avNone(i)));
    mp(i)=avLMAN(i)+avNone(i)/2;
end
chisqLM=chisqLM/(ending-starting+1); %normalize by number of points

chisqLM=chisqLM/mean(mp); %normalize by mean pitch