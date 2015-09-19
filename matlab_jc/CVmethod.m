function [gpre,gpost]=CVmethod(Data,Datakeep,index,width)
% Method 4: use CV of individual notes to determine if LMAN is broadening
% or reducing variation
% e.g. CVmethod(BFLesion,[1 2 3 4 9],[950 450 550 450 300],25);
ptstart=index-(width*8)/2+1;
ptend=index+(width*8)/2;

for k=1:length(Datakeep)
    i=Datakeep(k);
    gpost(k)=(median(std(residuals(Data(i).UDpost20_1024(ptstart(k):ptend(k),:)))./abs(mean(residuals(Data(i).UDpost20_1024(ptstart(k):ptend(k),:))))));
    gpre(k)=(median(std(residuals(Data(i).UDpre20_1024(ptstart(k):ptend(k),:)))./abs(mean(residuals(Data(i).UDpre20_1024(ptstart(k):ptend(k),:))))));
end