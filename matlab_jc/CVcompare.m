function CVs=CVcompare(p05p,f05p,p2post,f2post,ppre,fpre,p25pre,f25pre,first,second)

for i=0:5
    dist=2^(i+2)-4;

    a=first-dist;
    b=second+dist;


CV2post1=jc_plotabunch725(p05p,f05p,1000,a,b);
CV2post2=jc_plotabunch725(p2post,f2post,1000,a,b);
CV2pre1=jc_plotabunch725(ppre,fpre,1000,a,b);
CV2pre2=jc_plotabunch725(p25pre,f25pre,1000,a,b);

j=i+1;
CVs(j,1)=std(CV2pre1)/mean(CV2pre1);
CVs(j,2)=std(CV2pre2)/mean(CV2pre2);
CVs(j,3)=std(CV2post1)/mean(CV2post1);
CVs(j,4)=std(CV2post2)/mean(CV2post2);
end





