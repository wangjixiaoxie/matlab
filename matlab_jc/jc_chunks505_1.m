function jc_chunks505_1(baseline,shifted,inactivated,on,off,quantum)
quanta=(off-on)/10;
for i=1:length(baseline)
    for j=1:quanta
        beginning=(j-1)*(quantum)+1;
        ending=j*quantum;
        a(i,j)=median(baseline(i).pitches(beginning:ending));
    end
end
for j=1:quanta
    b(j)=mean(a(:,j));
end
pitcherbase=b;
%figure; hold on
plot(pitcherbase,'r'); 

%------------------------------------------------------------------
for i=1:length(shifted)
    for j=1:quanta
        beginning=(j-1)*(quantum)+1;
        ending=j*quantum;
        c(i,j)=median(shifted(i).pitches(beginning:ending));
    end
end
for j=1:quanta
    d(j)=mean(c(:,j));
end
pitchershift=d;
plot(pitchershift,'b'); 

%--------------------------------------------------------
for i=1:length(inactivated)
    for j=1:quanta
        beginning=(j-1)*(quantum)+1;
        ending=j*quantum;
        e(i,j)=median(inactivated(i).pitches(beginning:ending));
    end
end
for j=1:quanta
    f(j)=mean(e(:,j));
end
pitcherINACT=f;
plot(pitcherINACT,'k'); 