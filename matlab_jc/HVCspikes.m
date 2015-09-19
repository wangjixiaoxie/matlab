HVCspikes
% Channels (synchronized in time)
    % obs0r - song
    % obs1r - artificial FB
    % obs0 - could be a neural channel
    % obs1 - could be a neural channel
% Necessary to sort raw voltage profile into individual spikes
    % use SAUCY - find code in ~jsakata/matlab
    % stepsforSAUCY.m
% Take a look at song, spikes (see channels above)
    song=evsoundin('','o85pu54_vol7_2-070804.001.cbin','obs0r');
    spikes=evsoundin('','bl82bl81_vol7-013105.002.cbin','obs0');

    
    
    
    
    
%%% JC analysis %%%
    % Get raw data - goes 50ms into the past
    fvBspikes=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',1);
    fvBsong=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0r',1);

    % Calculate pitch contour
    for i=1:length(fvBsong) 
        shiftedBsong(i,:)=fvBsong(i).datt;
    end
    pitchBsong=jc_pitchmat1024(shiftedBsong,1024,1020,1,1500,2600,[1],'obs0',1);

    % Calculate rectified voltage contour and align to pitch contour
        for i=1:length(fvBspikes)
            shiftedBspikes(i,:)=fvBspikes(i).datt;
        end
        alignBspikes=resample(shiftedBspikes',1,4);
        alignBspikes=alignBspikes(128:end-128,:);
        rectBspikes=abs(alignBspikes);

    % Choose a fairly flat section of the pitch contours
        window=[420 520]; 

    % What is our premotor window in ms?
        premotorWindow=[-50:0];
    % The notes labeled as 'b' fall into two general categories
%         indb1=find(pitchbsong(256,:)>2060);
%         indb2=find(pitchbsong(256,:)<2060);
%     
    % Calculate correlations
    clear precorr 
    clear pretime
    for i=1:length(premotorWindow)
        pc=corrcoef(mean(rectBspikes(window+premotorWindow(i)*8,indb2)),median(pitchbsong(window,indb2)));
        precorr(i)=pc(2);
        pretime(i)=premotorWindow(i);
    end
    
    
%%%%  o85
    fvBspikes=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',1);
    fvBsong=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0r',1);

    % Calculate pitch contour
    for i=1:length(fvBsong) 
        shiftedBsong(i,:)=fvBsong(i).datt;
    end
    pitchBsong=jc_pitchmat1024(shiftedBsong,1024,1020,1,1500,2600,[1],'obs0',1);

    % Calculate rectified voltage contour and align to pitch contour
        for i=1:length(fvBspikes)
            shiftedBspikes(i,:)=fvBspikes(i).datt;
        end
        alignBspikes=resample(shiftedBspikes',1,4);
        alignBspikes=alignBspikes(128:end-128,:);
        rectBspikes=abs(alignBspikes);


ind1=find(pitchBsong(450,:)<2000);
clear spikings count
for i=1:size(alignBspikes,2)
spikings(i).data=find(alignBspikes(:,i)<(mean(alignBspikes(:,i))-2*std(alignBspikes(:,i))));
count(i)=length(find(spikings(i).data>80 & spikings(i).data<300));
end  
corrcoef(count(ind1),mean(pitchBsong(470:520,ind1)))
%%%%

    fvMspikes=findwnoteJC('batchnotes','m','','',0,[6000 8100],8500,1,'obs0',1);
    fvMsong=findwnoteJC('batchnotes','m','','',0,[6000 8100],8500,1,'obs0r',1);

    % Calculate pitch contour
    for i=1:length(fvMsong) 
        shiftedMsong(i,:)=fvMsong(i).datt;
    end
    pitchMsong=jc_pitchmat1024(shiftedMsong,1024,1020,1,1500,2600,[1],'obs0',1);

        for i=1:length(fvMspikes)
            shiftedMspikes(i,:)=fvMspikes(i).datt;
        end
        alignMspikes=resample(shiftedMspikes',1,4);
        alignMspikes=alignMspikes(128:end-128,:);
        rectMspikes=abs(alignMspikes);




clear spikings count
for i=1:size(alignMspikes,2)
spikings(i).data=find(alignMspikes(:,i)<(mean(alignMspikes(:,i))-2*std(alignMspikes(:,i))));
count(i)=length(find(spikings(i).data>1 & spikings(i).data<200));
end  
corrcoef(count,mean(pitchMsong(340:400,:)))

%%%%%%%
%%%%% pu77bk41
    fvBspikes=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',1);
    fvBsong=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0r',1);

    % Calculate pitch contour
    for i=1:length(fvBsong) 
        shiftedBsong(i,:)=fvBsong(i).datt;
    end
    pitchBsong=jc_pitchmat1024(shiftedBsong,1024,1020,1,1500,2600,[1],'obs0',1);

    % Calculate rectified voltage contour and align to pitch contour
        for i=1:length(fvBspikes)
            shiftedBspikes(i,:)=fvBspikes(i).datt;
        end
        alignBspikes=resample(shiftedBspikes',1,4);
        alignBspikes=alignBspikes(128:end-128,:);
        rectBspikes=abs(alignBspikes);
clear spikings count
for i=1:size(alignBspikes,2)
spikings(i).data=find(alignBspikes(:,i)>(mean(alignBspikes(:,i))+1*std(alignBspikes(:,i))));
count(i)=length(find(spikings(i).data>1 & spikings(i).data<160));
end  
corrcoef(count,mean(pitchBsong(240:400,:)))

