load /cardinal3/SyntaxBirds/FFsyntax1107.mat
%%%% r39g39 - FF - hamming window
datas=FFsyntax(3).data;
lenham=11;
h=hamming(lenham);
h=h*(lenham/sum(h));
for i=1:length(datas)
    acsfbase(i)=(datas(i).acsf & datas(i).baseline);
    acsfwn(i)=(datas(i).acsf & ~datas(i).baseline);
    apv(i)=(~datas(i).acsf);
end
acsfwnpitch=[];
acsfwntimes=[];
for i=1:length(datas)
    if acsfwn(i);
        times=timing3(datas(i).fvals);
        pts=mean(datas(i).pitch(datas(1).window,:));
        acsfwntimes=[acsfwntimes times];
        acsfwnpitch=[acsfwnpitch pts];
    end
end
figure;hold on;
unidays=find(diff(acsfwntimes)>8);
unistart=[1 unidays+1];
uniend=[unidays length(acsfwntimes)];
for i=1:length(unistart)
    times=acsfwntimes(unistart(i):uniend(i));
    pts=acsfwnpitch(unistart(i):uniend(i));
    hamfilttimes=[];hamfiltpts=[];
    for ii=1:length(times)-lenham
        hamfilttimes(ii)=mean(times(ii:ii+lenham-1).*h');
        hamfiltpts(ii)=mean(pts(ii:ii+lenham-1).*h');
    end
    plot(hamfilttimes,hamfiltpts,'k.','Markersize',15)
end
for i=find(acsfbase+apv)
    times=mean(timing3(datas(i).fvals));
    pts=mean(datas(i).pitch(datas(1).window,:));
    ptmn(i)=mean(mean(datas(i).pitch(datas(1).window,:)));
    guess=[];
    for ii=1:1000
        randsam=ceil(length(pts)*rand(1,length(pts)));
        guess(ii)=mean(pts(randsam));
    end
    pt05=prctile(guess,5);
    pt95=prctile(guess,95);
%     times=timing3(datas(i).fvals);
%     pts=mean(datas(i).pitch(datas(1).window,:));
%     hamfilttimes=[];hamfiltpts=[];
%     for ii=1:length(datas(i).fvals)-lenham
%         hamfilttimes(ii)=mean(times(ii:ii+lenham-1).*h');
%         hamfiltpts(ii)=mean(pts(ii:ii+lenham-1).*h');
%     end
    if acsfbase(i)
        plot(times,ptmn(i),'b.','Markersize',15)
        plot([times times],[pt05 pt95],'b-','Linewidth',2)
    else
        plot(times,ptmn(i),'r.','Markersize',15)
        plot([times times],[pt05 pt95],'r-','Linewidth',2)
    end
end
plot([6080 6200],[mean([ptmn(1) ptmn(3)]) mean([ptmn(1) ptmn(3)])],'k')



%%%% pu67bk2 - FF
datas=FFsyntax(1).data;
indAPVpoint=[2 4]
indACpoint=1;
indACcurve=3;
figure;hold on;
ravgwin=10;
for i=1:length(datas)
    if i==indACcurve
            [num den]=butter(5,1/20,'low');
            y=filtfilt(num,den,mean(datas(i).pitch(datas(1).window,:)));      
            plot(timing3(datas(i).fvals),y,'b')
    else
        plot(runningmedian(timing3(datas(i).fvals),ravgwin),runningaverage(mean(datas(i).pitch(datas(1).window,:)),ravgwin),'r*')
    end
end

%%%% r37g7 - FF
datas=FFsyntax(2).data;
indAPVpoint=[2 5 7]
indACpoint=[1 3];
indACcurve=[4 6];
figure;hold on;
ravgwin=10;
for i=1:length(datas)
    if ~isempty(find(indACcurve==i))
            [num den]=butter(5,1/20,'low');
            y=filtfilt(num,den,mean(datas(i).pitch(datas(1).window,:)));      
            plot(timing3(datas(i).fvals),y,'b')
    else
        %plot(runningmedian(timing3(datas(i).fvals),ravgwin),runningaverage(mean(datas(i).pitch(datas(1).window,:)),ravgwin),'r*')
    end
end    

