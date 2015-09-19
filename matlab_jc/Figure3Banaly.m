function [BSlopesLong,BSlopesShort]=Figure3Banaly(Alldata,prctiles)
% This does sensitivity analysis for Figure 3 - beginning/middle/end and
% at different contingencies.
% prctiles=[50 55 60 65 70 75 80 85 90];

for i=1:length(prctiles)
    for j=1:3 % begin,middle,end
        BaselineSlopes=Figure3analy(j,Alldata,prctiles(i));
        BSlopesLong(i,j)=mean(BaselineSlopes(Alldata(1).ind_longnotes));
        BSlopesShort(i,j)=mean(BaselineSlopes(Alldata(1).ind_shortnotes));
    end
end

g=8;


