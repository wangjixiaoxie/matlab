function jcsim(FFstats13pm)
for i=1:1000
a=rand(1,325);
a=round(a*648)+1;
aa=median(FFstats13pm(a))-median(FFstats13pm);
cc(i)=2*aa;
end
g=5;