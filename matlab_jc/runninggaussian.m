function [dataOUT]=runninggaussian(dataIN,windowsize)
gg=gausswin(windowsize);
for i=1:length(dataIN)-windowsize
    dataOUT(i:i+windowsize-1)=dataIN(i:i+windowsize-1)*gg;
end
dataOUT=dataOUT*(1/sum(gg));