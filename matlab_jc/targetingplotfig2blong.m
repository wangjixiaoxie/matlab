function targetingplotfig2blong(Alldata,notesizeshiftms)

% prctile_cutoff=0.8;
% tmp_cutoff=24;
% shiftedmethod=1;
% notesizeshiftms=0;
notesizeshift=notesizeshiftms*8;
firstover1sig=[8 9 9 5 13 13 11 4 7 7 9 7 6 11 12 8];


for i=1:length(Alldata)
    ons=Alldata(i).startnote;
    offs=Alldata(i).endnote;
    exp=Alldata(i).exp;
    shift_direction=Alldata(i).exp(1).direction;
    
    % get an indicator of the mean pitch in each 4-14hr group of song
    count=0;

    for j=firstover1sig(i)
            z(i).toffs=exp(j).toffset;
      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Alldata)
    ons(i)=Alldata(i).startnote+50;
    offs(i)=Alldata(i).endnote-50;
end


mm(1).dat=mean(Alldata(1).exp(8).selectedpitchcurves(ons(1):offs(1),1:20)')-mean(Alldata(1).baselineAC(ons(1):offs(1))');
mm(2).dat=mean(Alldata(3).exp(9).selectedpitchcurves(ons(3):offs(3),:)')-mean(Alldata(3).baselineAC(ons(3):offs(3))');
mm(3).dat=mean(Alldata(7).exp(11).selectedpitchcurves(ons(7):offs(7),:)')-mean(Alldata(1).baselineAC(ons(7):offs(7))');
mm(4).dat=mean(Alldata(8).exp(4).selectedpitchcurves(ons(8):offs(8),:)')-mean(Alldata(1).baselineAC(ons(8):offs(8))');
mm(5).dat=mean(Alldata(9).exp(7).selectedpitchcurves(ons(9):offs(9),:)')-mean(Alldata(1).baselineAC(ons(9):offs(9))');
mm(6).dat=mean(Alldata(12).exp(7).selectedpitchcurves(ons(12):offs(12),1:25)')-mean(Alldata(1).baselineAC(ons(12):offs(12))');
mm(7).dat=mean(Alldata(13).exp(6).selectedpitchcurves(ons(13):offs(13),:)')-mean(Alldata(1).baselineAC(ons(13):offs(13))');
mm(8).dat=mean(Alldata(14).exp(11).selectedpitchcurves(ons(14):offs(14),:)')-mean(Alldata(1).baselineAC(ons(14):offs(14))');
mm(9).dat=mean(Alldata(15).exp(12).selectedpitchcurves(ons(15):offs(15),:)')-mean(Alldata(1).baselineAC(ons(15):offs(15))');
mm(10).dat=mean(Alldata(16).exp(8).selectedpitchcurves(ons(16):offs(16),1:40)')-mean(Alldata(1).baselineAC(ons(16):offs(16))');
for i=1:10
    mm(i).dat=abs(mm(i).dat/max(abs(mm(i).dat)));
    [maxi(i),indi(i)]=max(mm(i).dat);
end

abc=zeros(length(Alldata),1000);

for n=1:length(mm);   % JUST THE LONG NOTES
left=indi(n)-1;
right=length(mm(n).dat)-indi(n);
abc(n,1)=mm(n).dat(indi(n));
if right>left
    for m=2:left
        abc(n,m)=mean([mm(n).dat(indi(n)-m) mm(n).dat(indi(n)+m)]);
    end
    for m=left+1:right
        abc(n,m)=mm(n).dat(indi(n)+m);
    end
else
    for m=2:right
        abc(n,m)=mean([mm(n).dat(indi(n)-m) mm(n).dat(indi(n)+m)]);
    end
    for m=right+1:left
        abc(n,m)=mm(n).dat(indi(n)-m);
    end
end
end
    
for nn=1:size(abc,2)
    indices=find(abc(:,nn)~=0);
    avcurve(nn)=mean(abc(indices,nn));
    avup(nn)=mean(abc(indices,nn))+std(abc(indices,nn));
    avdn(nn)=mean(abc(indices,nn))-std(abc(indices,nn));
    xaxis(nn)=nn/8;
end
plot(xaxis,avcurve);hold on;plot(xaxis,avup);hold on;plot(xaxis,avdn)
xlim([0 85])

% Get all toffsets and graph on the same plot
toall=[];
for ii=1:length(Alldata(1).ind_longnotes)
    i=Alldata(1).ind_longnotes(ii);
    toall=[toall;z(i).toffs-mean(z(i).toffs)];
end
histo=hist(abs(toall)/8,50);
%hold on;plot(histo/max(histo),'r')

% Do it for right over 1sigma


