for i=51:200
    resid=pitch50r(700:1000,i)-mean(pitch50r(700:1000,1:200)')';
    indesc=find(allescapes50r(1:i-1)==1);
    indhit=find(allescapes50r(1:i-1)==0);
    lastescs=max(indesc);
    lasthit=max(indhit);
    desc(i)=sum(abs(pitch50r(700:1000,i)-(pitch50r(700:1000,indesc(end-1)))))/301;
    dhit(i)=sum(abs(pitch50r(700:1000,i)-(pitch50r(700:1000,indhit(end-1)))))/301;
    dres(i)=sum(abs(resid))/301;
end