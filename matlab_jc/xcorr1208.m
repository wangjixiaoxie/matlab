function [slopesAC,slopesINA]=xcorr1208(Alldata2)
g=[1 4 5 6 10 11 15 16];
for i=1:length(g)
    clear xcAC
    clear xcINA
    clear pitchAC
    clear pitchINA
    gg=g(i);
    pitchAC1=Alldata2(gg).baselineAC(Alldata2(gg).basevaronset:Alldata2(gg).basevaroffset-40,1:77);
    pitchINA1=Alldata2(gg).baselineINA(Alldata2(gg).basevaronset:Alldata2(gg).basevaroffset-40,1:77);
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
    [peakINA,index]=max(median(xcINA'));
    [peakAC,index]=max(median(xcAC'));
    xcINAn=median(xcINA')./peakINA;
    xcACn=median(xcAC')./peakAC;
    
    hold on;plot(median(xcAC'))
    hold on;plot(median(xcINA'),'r')
   peaksINA(i)=peakINA;
   peaksAC(i)=peakAC;
   storeINA(i,:)=xcINAn(index-100:index+100);
   storeAC(i,:)=xcACn(index-100:index+100);
end

%%%% Timescale (slope of variation)
    %%%% Get the slopes
    x=(0:1/8:39/8);
    for i=1:size(storeAC,1)
        kk1=polyfit(x,storeAC(i,101:140),1);
        kk2=polyfit(x,storeINA(i,101:140),1);
        slopesAC(i)=-1*kk1(1);
        slopesINA(i)=-1*kk2(1);
    end

    %%%% Plot the raw data
   figure;plot(1,slopesAC,'*','Color','k');hold on;plot(2,slopesINA,'*','Color','k')
    hold on;plot(1,slopesAC([2 3]),'*','Color','r')
    hold on;plot(2,slopesINA([2 3]),'*','Color','r')
     for i=1:length(g)
         hold on;plot([1 2], [slopesAC(i) slopesINA(i)], '-','Color','k')
     end


%%%% Magnitude
    figure;plot(1,peaksAC,'*','Color','k');hold on;plot(2,peaksINA,'*','Color','k')
    hold on;plot(1,peaksAC([2 3]),'*','Color','r')
    hold on;plot(2,peaksINA([2 3]),'*','Color','r')
     for i=1:length(g)
         hold on;plot([1 2], [peaksAC(i) peaksINA(i)], '-','Color','k')
     end
