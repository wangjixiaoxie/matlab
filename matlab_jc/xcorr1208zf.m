function xcorr1208zf(AlldataZFlesion)
g=[1 2 3 4 5 6 7];
for i=1:6
    clear xcAC
    clear xcINA
    clear pitchAC
    clear pitchINA
    gg=g(i);
    pitchAC1=AlldataZFlesion(gg).pitchUDpre(AlldataZFlesion(gg).basevaronset+40:AlldataZFlesion(gg).basevaroffset-40,1:25);
    pitchINA1=AlldataZFlesion(gg).pitchUDpost(AlldataZFlesion(gg).basevaronset+40:AlldataZFlesion(gg).basevaroffset-40,1:25);
    for jj=1:size(pitchAC1,2)
        pitchAC(:,jj)=pitchAC1(:,jj)-mean(pitchAC1')';
    end
    for kk=1:size(pitchINA1,2)
        pitchINA(:,kk)=pitchINA1(:,kk)-mean(pitchINA1')';
    end
    for j=1:size(pitchAC,2)
        xcAC(:,j)=xcorr(pitchAC(:,j));
    end
    for k=1:size(pitchINA,2)
        xcINA(:,k)=xcorr(pitchINA(:,k));
    end
    hold on;plot(median(xcAC')./max(median(xcAC')))
    hold on;plot(median(xcINA')./max(median(xcINA')),'r')
end