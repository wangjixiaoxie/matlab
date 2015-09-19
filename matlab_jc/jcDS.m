%for i=1:length(DShifts)
    i=5;
    clear dayvals;
    clear indhit
   % day1=DShifts(i).fvALL(1).month*30+DShifts(i).fvALL(1).day;
    lengthbase=size(DShifts(i).pitchBaseline,2);
    for ii=1:length(DShifts(i).fvALL)
        dayvals(ii)=DShifts(i).fvALL(ii).month*30+DShifts(i).fvALL(ii).day;
        indhit(ii)=DShifts(i).fvALL(i).TRIG;
        if isequal(DShifts(i).fvALL(ii).lbl(DShifts(i).fvALL(ii).ind),'b')
            indhit(ii)=0;
        else
            indhit(ii)=1;
        end
                %indhit(ii)=DShifts(i).fvALL(ii).TRIG;
    end
    dayvals=dayvals(lengthbase+1:end);
    indhit=indhit(lengthbase+1:end);
    minD=min(dayvals);
    maxD=max(dayvals);
    pitchAEsc=[];
    pitchBEsc=[];
    avEsc=[];
    avHit=[];
    avAll=[];
    pitchAHit=[];
    pitchBHit=[];
    
    for j=minD:maxD
        indDay=find(dayvals==j);
        avAll(j,:)=mean(DShifts(i).pitchALL(:,indDay+lengthbase)');
        indEday=find(indhit(indDay)==0);
        indHday=find(indhit(indDay)==1);
        avEsc(j,:)=mean(DShifts(i).pitchALL(:,indEday+lengthbase)');
        avHit(j,:)=mean(DShifts(i).pitchALL(:,indHday+lengthbase)');
        medtoffs=median(DShifts(i).toffset);
        pitchAEsc(j)=median(avEsc(j,medtoffs-64:medtoffs));
        pitchBEsc(j)=median(avEsc(j,medtoffs-64+192:medtoffs+192));
        pitchAHit(j)=median(avHit(j,medtoffs-64:medtoffs));
        pitchBHit(j)=median(avHit(j,medtoffs-64+192:medtoffs+192));
        pitchAavg(j)=median(avAll(j,medtoffs-64:medtoffs));
        pitchBavg(j)=median(avAll(j,medtoffs-64+192:medtoffs+192));
    end
    if isequal(DShifts(i).dirA,'up')
        dirA=1;
        dirB=-1;
    else
        dirA=-1;
        dirB=1;
    end
figure;plot(dirA*(pitchAEsc(minD:maxD-1)-pitchAHit(minD:maxD-1)),dirA*(pitchAavg(minD+1:maxD)-pitchAavg(minD:maxD-1)),'*')
hold on;plot(dirB*(pitchBEsc(minD:maxD-1)-pitchBHit(minD:maxD-1)),dirB*(pitchBavg(minD+1:maxD)-pitchBavg(minD:maxD-1)),'*','Color','r')

j=7;
diff=[];
[b,ix]=sort(timing3(DShifts(j).fvALL));
for i=1:size(DShifts(j).pitchALL,2)
    diff(:,i)=DShifts(j).pitchALL(DShifts(j).onset:DShifts(j).offset,i)'-mean(DShifts(j).pitchBaseline(DShifts(j).onset:DShifts(j).offset,:)');
end
figure;image(diff(:,ix))
