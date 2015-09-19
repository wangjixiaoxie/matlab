function pitchT=jc_exclude(pitchA)
count=1;
init=220;
fin=455;
width=fin-init+1;
for i=1:size(pitchA,2)
    for j=init:fin
        Adt(j)=abs(pitchA(j-1,i)-pitchA(j,i));
    end
    if max(Adt)<20 %&& min(normA(init:fin,i))>-0.1% &&max(normA(init:fin,i))<0.06
        pitchT(:,i)=pitchA(:,i);
        count=count+1;
    end
end
