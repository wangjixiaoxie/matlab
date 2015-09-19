function aligned=TraceAlign(smooth)
avgnote=mean(smooth);
for i=1:size(smooth,1)
    h=xcov(avgnote(1:7000),smooth(i,1:7000));
    [peak,index]=max(h);
    shift(i)=(index-7000);
    starter=500-round(shift(i));
    ender=starter+6000;
    aligned(i,:)=smooth(i,starter:ender); 
end