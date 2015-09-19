figure
FFrecov=[260 450 1 50 50 50 0 0 0 50 0 50 140 1 1 40 1];
for i=[1:2 4:6 10 12 13:17]
    target=round(median(Experiment(i).on:Experiment(i).off));
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    [x,sorted]=sort(Experiment(i).timeACpre);
    t1=[1/length(sorted):1/length(sorted):1];
    baseline=median(Experiment(i).pitchACpre(target,end-50:end));
    if i==16
        baseline=median(Experiment(i).pitchACpre(target,:));
    end
    hold on;plot(runningaverage(t1,round(length(t1)/20)),coef*(runningaverage(Experiment(i).pitchACpre(target,sorted),round(length(t1)/20))-baseline),'-')
    ratime(i,1).data=runningaverage(t1,round(length(t1)/20));
    raFF(i,1).data=coef*(runningaverage(Experiment(i).pitchACpre(target,sorted),round(length(t1)/20))-baseline);
    [x,sorted2]=sort(Experiment(i).timeAPV);
    t2=[1+1/length(sorted2):1/length(sorted2):2];
    hold on;plot(runningaverage(t2,round(length(t2)/20)),coef*(runningaverage(Experiment(i).pitchAPV(target,sorted2),round(length(t2)/20))-baseline),'-','Color','g')
    ratime(i,2).data=runningaverage(t2,round(length(t2)/20));
    raFF(i,2).data=coef*(runningaverage(Experiment(i).pitchAPV(target,sorted2),round(length(t2)/20))-baseline);
    [x,sorted3]=sort(Experiment(i).timeAPVwn);
    t3=[2+1/length(sorted3):1/length(sorted3):3];
    hold on;plot(runningaverage(t3,round(length(t3)/20)),coef*(runningaverage(Experiment(i).pitchAPVwn(target,sorted3),round(length(t3)/20))-baseline),'-','Color','r')
    ratime(i,3).data=runningaverage(t3,round(length(t3)/20));
    raFF(i,3).data=coef*(runningaverage(Experiment(i).pitchAPVwn(target,sorted3),round(length(t3)/20))-baseline);
    [x,sorted4]=sort(Experiment(i).timeACpost(FFrecov(i):end));
    t4=[3+1/length(sorted4):1/length(sorted4):4];
    time4(i).data=x;
    hold on;plot(runningaverage(t4,round(length(t4)/20)),coef*(runningaverage(Experiment(i).pitchACpost(target,sorted4),round(length(t4)/20))-baseline),'-','Color','k')
    ratime(i,4).data=runningaverage(t4,round(length(t4)/20));
    raFF(i,4).data=coef*(runningaverage(Experiment(i).pitchACpost(target,sorted4),round(length(t4)/20))-baseline);
end


ind=[1:2 10  12 13:16];
for i=[1:2 10 12 13:16]
    for j=1:20
        width=length(raFF(i,1).data)/20;
        mnFF1(i,j)=mean(raFF(i,1).data(1+(j-1)*width:j*width));
        width=length(raFF(i,2).data)/20;
        mnFF2(i,j)=mean(raFF(i,2).data(1+(j-1)*width:j*width));
        width=length(raFF(i,4).data)/20;
        mnFF4(i,j)=mean(raFF(i,4).data(1+(j-1)*width:j*width));
        width=length(raFF(i,3).data)/20;
        mnFF3(i,j)=mean(raFF(i,3).data(1+(j-1)*width:j*width));
    end
end
figure;hold on;
plot([1:20],mean(mnFF1(ind,:)))
plot([21:40],mean(mnFF2(ind,:)),'g')
plot([41:60],mean(mnFF3(ind,:)),'r')
plot([61:80],mean(mnFF4(ind,:)),'k')



ind=[1:2 4:6 10 12:17]
for i=ind
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    targ=Experiment(i).on;
    baseline=mean(Experiment(i).pitchACpre(targ,find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-5)));
    [x,sortedAPV]=sort(Experiment(i).timeAPV);
    APVbegin(i)=coef*(mean(Experiment(i).pitchAPV(targ,sortedAPV(end-15:end)))-baseline);
    APVend(i)=coef*(mean(Experiment(i).pitchAPVwn(targ,end-15:end))-baseline);
    for j=1:70
        times1=find(Experiment(i).timeACpost>min(Experiment(i).timeACpost)+j-1);
        times2=find(Experiment(i).timeACpost(times1)<min(Experiment(i).timeACpost)+j);
        place=[];
        place=coef*(Experiment(i).pitchACpost(targ,times1(times2))-baseline);
        if ~isempty(place)
            ACpost(i,j)=mean(place);
        else
            ACpost(i,j)=0;
        end
    end
end
for i=1:70
    mnACpost(i)=mean(ACpost(find(ACpost(:,i)~=0),i));
    sdACpost(i)=std(ACpost(find(ACpost(:,i)~=0),i));
end

figure;plot(mnACpost+sdACpost/sqrt(12))
hold on;plot(mnACpost-sdACpost/sqrt(12))
hold on;plot([-10 0],[mean(APVbegin(ind));mean(APVend(ind))])
hold on;plot([-10 0],[mean(APVbegin(ind))+std(APVbegin(ind))/sqrt(12);mean(APVend(ind))+std(APVend(ind))/sqrt(12)])
hold on;plot([-10 0],[mean(APVbegin(ind))-std(APVbegin(ind))/sqrt(12);mean(APVend(ind))-std(APVend(ind))/sqrt(12)])




figure;plot([APVbegin;APVend])

figure;hold on
for i=[1:2 4:6 10 12:17]
if isequal(Experiment(i).DIR,'up');
coef=1;
else
coef=-1;
end
hold on;plot(runningaverage(Experiment(i).timeACpost,15)-min(Experiment(i).timeACpost),coef*(runningaverage(Experiment(i).pitchACpost(Experiment(i).on,:),15)-mean(Experiment(i).pitchACpre(Experiment(i).on,end-90:end))),'k')
end



%%%%
% with sleep
%%%%
ind=[1:2 4:6 10 12:17]
for i=ind
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    targ=Experiment(i).on;
    baseline=mean(Experiment(i).pitchACpre(targ,find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-5)));
    [x,sortedAPV]=sort(Experiment(i).timeAPV);
    APVbegin(i)=coef*(mean(Experiment(i).pitchAPV(targ,sortedAPV(end-15:end)))-baseline);
    APVend(i)=coef*(mean(Experiment(i).pitchAPVwn(targ,end-15:end))-baseline);
    % find first night of sleep
    [x,sortACpost]=sort(Experiment(i).timeACpost);
    clear gaps
    for k=2:length(Experiment(i).timeACpost)
        gaps(k-1)=Experiment(i).timeACpost(sortACpost(k))-Experiment(i).timeACpost(sortACpost(k-1));
    end
    nights=find(gaps>8);
    nighttime1=round((Experiment(i).timeACpost(sortACpost(nights(1)-1))-min(Experiment(i).timeACpost)));
    g(i)=nighttime1;
    for j=1:nighttime1
        times1=find(Experiment(i).timeACpost>min(Experiment(i).timeACpost)+j-1);
        times2=find(Experiment(i).timeACpost(times1)<min(Experiment(i).timeACpost)+j);
        place=[];
        place=coef*(Experiment(i).pitchACpost(targ,times1(times2))-baseline);
        if ~isempty(place)
            ACpost(i,j)=mean(place);
        else
            ACpost(i,j)=0;
        end
    end
    for j=nighttime1:70
        times1=find(Experiment(i).timeACpost>min(Experiment(i).timeACpost)+j-1);
        times2=find(Experiment(i).timeACpost(times1)<min(Experiment(i).timeACpost)+j);
        place=[];
        place=coef*(Experiment(i).pitchACpost(targ,times1(times2))-baseline);
        if ~isempty(place)
            ACpost2(i,j)=mean(place);
        else
            ACpost2(i,j)=0;
        end

    end
    
end
for i=1:70
    mnACpost2(i)=mean(ACpost2(find(ACpost2(:,i)~=0),i));
    sdACpost2(i)=std(ACpost2(find(ACpost2(:,i)~=0),i));
end

figure;plot(mnACpost+sdACpost/sqrt(12))
hold on;plot(mnACpost-sdACpost/sqrt(12))
hold on;plot([-10 0],[mean(APVbegin(ind));mean(APVend(ind))])
hold on;plot([-10 0],[mean(APVbegin(ind))+std(APVbegin(ind))/sqrt(12);mean(APVend(ind))+std(APVend(ind))/sqrt(12)])
hold on;plot([-10 0],[mean(APVbegin(ind))-std(APVbegin(ind))/sqrt(12);mean(APVend(ind))-std(APVend(ind))/sqrt(12)])

%%%%%%%%%%
%%%%%%%%%%
%% 12.9.09
%%%%%%%%%%
% Figure 1 - Adaptive shifts (entire next day post - pvs 50 songs pre)
figure;plot(1,shiftvals,'v','Color','b','MarkerSize',10)
hold on;plot([0 2],[0 0],'k')
hold on;plot([0.9 1.1],[mean(shiftvals) mean(shiftvals)],'k')
xlim([0.8 1.2])
hold on;plot(1,shiftvals([2 4 7 10 13]),'v','Color','r','MarkerSize',10)
% Figure 1prime - Adaptive shifts 2 
%Same as above but with pre being all data in the pre file - better and just as significant
shiftvals2=shiftvals;
shiftvals2(1:10)=shiftvals2(1:10)-differ([1:2 4:6 10 12:13 15:16]);
figure;plot(1,shiftvals2,'v','Color','b','MarkerSize',10)
hold on;plot([0 2],[0 0],'k')
hold on;plot([0.9 1.1],[mean(shiftvals2) mean(shiftvals2)],'k')
xlim([0.8 1.2])
hold on;plot(1,shiftvals2([2 4 7 10 13]),'v','Color','r','MarkerSize',10)
%%% What is the change with AP5 on? - subtract end of AP5 post (last 20%) - from end of AP5 w/o WN (last 20%)
    for i=[1:2 4:6 10 12:13 15:16]
    AP5n(i)=mean(raFF(i,3).data(end-length(raFF(i,3).data)/5:end))-mean(raFF(i,2).data(end-length(raFF(i,2).data)/5:end));
    end
    AP5n(17)=mean(pitch1201apvwn(240,end-size(pitch1201apvwn,2)/5:end))-mean(pitch1201apv(240,end-size(pitch1201apv,2)/5:end));
    AP5n(18)=mean(pitch1203apvwnC(260,end-size(pitch1203apvwnC,2)/5:end))-mean(pitch1203apvC(260,end-size(pitch1203apvC,2)/5:end));
    AP5n(19)=-1*(mean(pitch1206apvwnC(260,end-size(pitch1206apvwnC,2)/5:end))-mean(pitch1206apvC(260,end-size(pitch1206apvC,2)/5:end)));
    AP5n=AP5n([1:2 4:6 10 12:13 15:19]);
hold on;plot(0.5,AP5n,'v','Color','b','MarkerSize',10)
xlim([0.3 1.2])
plot([0.4 0.6],[mean(AP5n) mean(AP5n)],'k')
plot(0.5,AP5n([2 4 7 10 13]),'v','Color','r','MarkerSize',10)
plot([0.5;1],[AP5n;shiftvals2],'k')


%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%% 12.15.09 %%%

%%% Control notes - AP5 and WN but no learning signal %%%
% Experiments 1:2 10 12:13 16 use the next note in the repeat as the ctl.
% Experiment 15 uses the long low stack (whereas the short high stack is
    % targeted).
% Experiments 18:20 use the low stack or high stack (whichever one is not
    % targeted).
ind=[1:2 10 12:13 15 16 18:20]
for i=ind
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    targ=Experiment(i).onCTL;
    baseline=mean(Experiment(i).pitchACpreCTL(targ,find(Experiment(i).timeACpreCTL>Experiment(i).timeACpreCTL(end)-5)));
    % find first night of sleep
        [x,sortACpost]=sort(Experiment(i).timeACpostCTL);
        clear gaps
        gaps(1)=Experiment(i).timeACpostCTL(sortACpost(1))-Experiment(i).timeAPVwn(end);
        for k=2:length(Experiment(i).timeACpostCTL)
            gaps(k)=Experiment(i).timeACpostCTL(sortACpost(k))-Experiment(i).timeACpostCTL(sortACpost(k-1));
        end
        nights=find(gaps>8);
        if length(nights)==1
            indexed=sortACpost(nights(1):end);
        else
            indexed=sortACpost(nights(1):nights(2));
        end
        diffCTL(i)=coef*(mean(Experiment(i).pitchACpostCTL(targ,indexed))-baseline);
end
for i=1:70
    mnACpost2(i)=mean(ACpost2(find(ACpost2(:,i)~=0),i));
    sdACpost2(i)=std(ACpost2(find(ACpost2(:,i)~=0),i));
end
