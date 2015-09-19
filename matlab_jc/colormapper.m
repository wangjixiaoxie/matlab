function colormapper(Alldata,targPositions)
    dataMatrix=zeros(5,100);
for i=2:length(Alldata)
    dat=Alldata(i).exp.SDtrace;
    dat=dat(Alldata(i).exp.notonset:Alldata(i).exp.notoffset);
    if Alldata(i).exp.max==1
        a=max(dat);
    else
        a=min(dat);
    end

    dat=dat/a;
    dat=dat*64; %normalize based on %of extremum
    dat=resample(dat,100,length(dat));

    histo=ceil((targPositions(i)-0.2)*5+1);
    if sum(dataMatrix(histo,:))==0
        dataMatrix(histo,:)=dat;
    else 
        dataMatrix(histo,:)=(dataMatrix(histo,:)+dat)./2;
    end
end
figure;image(dataMatrix);colormap(hot)