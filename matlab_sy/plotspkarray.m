function plotspkarray(spkarray,hght,offset)
        spks2=spkarray'
        spks2(:,2)=hght;
        plotrasters(spks2,hght,offset);