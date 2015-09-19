%try to relate the start_time of the response to the last syllable
%transition_time


    clear array;
    mnrsptm=mean(conscresp,2);
    mnrsptm=mnrsptm/(1/binsize)
    
    for i=1:length(mnrsptm)
    %plot rasters for syllable spacing
        array{i}(:,1)=mnrsptm(i)-xlist;
        array{i}(:,2)=i;
        plotrasters3(array);
    end