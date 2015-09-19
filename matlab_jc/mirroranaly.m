% EvTaf version 3.1
    % Make sure it doesn't trigger on playbacks - it doesn't appear to do so.
    % Stims must be in a separate folder or it will crash.

% Get raw data
    fvals=findwnoteJC('batchJC2notes','a','','',0,[6500 8100],8500,1,'obs0',1);
    for i=1:length(fvals)
        shifted(i,:)=fvals(i).datt;
    end
    pitch=jc_pitchmat1024(shifted,1024,1020,2,2000,2700,[1],'obs0',1);
    
% Set parameters    
    wnplacement=850; % Choose based on pitch contours
    delaywnplacement=(wnplacement+128)/8; % 128pts=512sampling points for the window half-width of 1024 pts
    cutoff=prctile(mean(pitch(wnplacement-64:wnplacement,:)),30); % "artificial escapes"
    indwn=find(mean(pitch(wnplacement-64:wnplacement,:))>cutoff); % syllables that artificially escape
    miscal=44; % -30ms (findwnoteJC) +10ms (miscalibration?)
% Make artificial WN
    % NOTE - use lowpasswn.m to make low pass wn for maximal verisimilitude
    % MAKE sure WN is at 32000Hz, not 44100Hz
    timewn=50; % in ms
    lengthwn=32000*(timewn/1000);
    rawdata=rand(lengthwn,1)*2-1;
    factorlouder=4; % How much WN is louder than syllable it's played within

% Make data for the stimuli songs
    for i=1:length(fvals)
        % For each new song, get the oscillogram
        if i==1
            count=1;
            [newdat(count).data,fs]=evsoundin('',fvals(i).fn,'obs0');
        else
            if ~isequal(fvals(i-1).fn,fvals(i).fn)
                count=count+1;
                [newdat(count).data,fs]=evsoundin('',fvals(i).fn,'obs0');
            end
        end
            onnumber=fvals(i).ind;
            noteonset=fvals(i).ons(fvals(i).ind)-miscal;
            rawdatascaled=rawdata*factorlouder*max(newdat(count).data(noteonset*32:noteonset*32+8500));
            zerobefore=(noteonset+delaywnplacement)*32; 
            zeroafter=size(newdat(count).data)-(zerobefore+lengthwn); 
            %addwn=[zeros(zerobefore,1);rawdatascaled;zeros(zeroafter,1)];
            if ~isempty(find(indwn==i))
                newdat(count).data(zerobefore:zerobefore+lengthwn-1)=newdat(count).data(zerobefore:zerobefore+lengthwn-1)+rawdatascaled;
            end
         % RUN THIS in the for loop to get the pitch data - make sure WN follows directly after the contingency
            % 256pts for the 8ms failure
            
               dtest=newdat(count).data(noteonset*32-256+64:noteonset*32-256+64+8500)';
               pitch2(:,i)=jc_pitchmat1024b(dtest,1024,1020,2,2000,2700,1,'obs0');
    end
    figure;plot(pitch,'k')
    hold on;plot(pitch2,'r')

% Open new folder for stim files
    for i=1:29
        name=['testdataug' num2str(i)];
        wavwrite(newdat(i).data/max(newdat(i).data),32000,16,name);
    end

    figure;plot(pitch,'k')
    hold on;plot(pitch2,'r')
    