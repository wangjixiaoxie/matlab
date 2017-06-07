function switches = makeswitchstructure(ng,transtime,switchtonotegroup,day,time)
%thistranstime acts as sign which switches to exclude

%day and time are datenums counting from 0, time is in ms

nswitches = sum(diff(ng)~=0);
switchix = find(diff(ng)~=0);
startcount = 1;

for i = 1:nswitches-1
    switches.preix{i} = startcount:switchix(i);
    switches.postix{i} = switchix(i)+1:switchix(i+1);
    switches.day(i) = day(switches.preix{i}(end));
    
    %get rid of renditions from previous day if notegroup stayed the same
    if day(switches.preix{i}(1))~=day(switches.preix{i}(end))
        fprintf('dropping block beginning on previous day; day %d \n',day(switches.preix{i}(1)))
    end
    if day(switches.postix{i}(end))~=day(switches.postix{i}(1))
%         fprintf('dropping block end on next day; first day was %d \n',day(switches.postix{i}(1)))
    end
    switches.preix{i}(day(switches.preix{i})~=day(switches.preix{i}(end))) = [];
    switches.postix{i}(day(switches.postix{i})~=day(switches.postix{i}(1))) = [];
    
    switches.preng(i) = mean(ng(switches.preix{i}));
    switches.postng(i) = mean(ng(switches.postix{i}));
    lastpretime = time(switches.preix{i}(end))./1000./60;
    firstposttime = time(switches.postix{i}(1))./1000./60;
    %damn, transtime ist nur auf die minute genau --> pretime - 60 um eine
    %minute abzurunden muesste reichen?
    %update fuer song switches: rechne 1ms auf file start time drauf falls
    %switch gleichzeitig quasi
    possibletranstimes = find(transtime.*60<=(firstposttime+(1/60/1000)) & transtime.*60>lastpretime-1);
    if length(possibletranstimes)==1 && switchtonotegroup(possibletranstimes)==switches.postng(i)
        switches.thistranstime(i) = transtime(possibletranstimes).*60;
    elseif isempty(possibletranstimes)
        fprintf('ACHTUNG: no transition time found %d ***************************************** \n',i) %big problem
        switches.thistranstime(i) = nan;

    elseif length(possibletranstimes)>1 && day(switches.preix{i}(end))~=day(switches.postix{i}(end))
        switches.thistranstime(i) = nan;
        fprintf('ignoring day transition notegorup %d to %d, day %d \n',switches.preng(i),switches.postng(i),day(switches.preix{i}(end)))    
    elseif length(possibletranstimes)==2
        assert(switchtonotegroup(possibletranstimes(end))==switches.postng(i))
        %middle block was left out - no song
        if switches.preng(i)==switchtonotegroup(possibletranstimes(end))
            disp('switching back to first NG') %weird situation but def no real switch
            switches.thistranstime(i) = nan;
        elseif (switches.preng(i)<2&&switchtonotegroup(possibletranstimes(1))<2) ||  (switches.preng(i)>=2&&switchtonotegroup(possibletranstimes(1))>=2)
            %left out block had same color as previous block
            %take second transition time
%             fprintf('left out block before switch, NG %d \n',switchtonotegroup(possibletranstimes(1)))
            switches.thistranstime(i) = transtime(possibletranstimes(end)).*60;
        elseif (switches.postng(i)<2&&switchtonotegroup(possibletranstimes(1))<2) ||  (switches.postng(i)>=2&&switchtonotegroup(possibletranstimes(1))>=2)
            %middle block had same color as right block
            %take first transition time
            switches.thistranstime(i) = transtime(possibletranstimes(1)).*60;
%             fprintf('left out block after switch, NG %d \n',switchtonotegroup(possibletranstimes(1)))
        else
            %change color twice: 0 2 1 oder 0 3 1
            switches.thistranstime(i) = nan;
            fprintf('changed light twice? %d %d %d \n',switches.preng(i),switchtonotegroup(possibletranstimes))
            %look into this if happens a lot
        end
            
    else
        disp('more than two transtime --- left out block -- figure this out')
        %usually happens mostly for fast switches when I start a new
        %experiment / test light
        switches.thistranstime(i)=nan;
    end
    
    startcount = switchix(i)+1;
end
%last switch
switches.preix{i+1} = startcount:switchix(i+1);
switches.postix{i+1} = switchix(i+1)+1:length(ng);
switches.day(i) = day(switches.preix{i+1}(end));

switches.preng(i+1) = mean(ng(switches.preix{i+1}));
switches.postng(i+1) = mean(ng(switches.postix{i+1}));
lastpretime = time(switches.preix{i+1}(end))./1000./60;
firstposttime = time(switches.postix{i+1}(1))./1000./60;
possibletranstimes = find(transtime.*60<firstposttime & transtime.*60>lastpretime-1);
if length(possibletranstimes)==1 && switchtonotegroup(possibletranstimes)==switches.postng(i+1)
    switches.thistranstime(i+1) = transtime(possibletranstimes).*60;
elseif isempty(possibletranstimes)
    disp('ACHTUNG: no transition time found')
    switches.thistranstime(i+1) = nan;
elseif length(possibletranstimes)>1 && day(switches.preix{i+1}(end))~=day(switches.postix{i+1}(end))
    switches.thistranstime(i+1) = nan;
    fprintf('ignoring day transition notegorup %d to %d \n',switches.preng(i),switches.postng(i))
else
    disp('more than one transtime --- left out block -- figure this out')
    switches.thistranstime(i+1) = nan;
end

fprintf('\n Number of switches %d, number of switches to analyze %d \n',length(switches.thistranstime),sum(~isnan(switches.thistranstime)))