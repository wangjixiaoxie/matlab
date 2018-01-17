function cleandirAuto(batch,wind,numwind,numnote)
% cleandirAuto(batch,wind,numwind,numnote,CHANSPEC)
% 
% wind is the window size in MS
% numwind is the number of notes in a window
% numnote is the numer of times numwind has to be observed in 1 file
%
% Threshold is calculated by mbatchamphist() based on distributiuon of 
% amplitude values across files in batch.

if (~exist('TH'))
    TH=5;
end

if (~exist('wind'))
    wind=1000;
end

wind = wind*1e-3; % convert to ms

if (~exist('numwind'))
    numwind=6;
end

if (~exist('numnote'))
    numnote=4;
end

% 
% if (~exist('CHANSPEC'))
%     CHANSPEC='obs0';
% end

% if (exist([batch '.keep']));
%     ans=input('.keep file exists, are you sure you want to proceed?  ','s');
%     if (ans=='y')
%     else
%         return;
%     end
% end

fid=fopen(batch,'r');
fkeep=fopen([batch,'.keep'],'w');
fdcrd=fopen([batch,'.dcrd'],'w');
disp(['working...']);

%calculate distribution of amplitudes in all songs in batch
[batchbins batchhist] = mbatchampdist(batch);
[pks,pksloc] = findpeaks(batchhist,'SORTSTR','descend');

threshold = 10^(batchbins(pksloc(1))+(batchbins(pksloc(2))-batchbins(pksloc(1))));

disp(['threshold = ' num2str(threshold)]);

while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end

   %disp(fn);

    [pth,nm,ext]=fileparts(fn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn,'0r');
        sm=evsmooth(dat,fs,0.01);
    elseif(strcmp(ext,'.cbin'))
        [dat,fs]=ReadCbinFile(fn, '0');
        sm=mquicksmooth(dat,fs);
    elseif(strcmp(ext,'.wav'))
        [dat,fs]=wavread(fn);
        sm=mquicksmooth(dat,fs);
    end
    %[ons,offs]=evsegment(sm,fs,5.0,30.0,TH);
    
    %threshold = mautothresh(fn,TH);
   
    [ons offs] = msegment(sm,fs,15,20,threshold);
    %filter vocalizations that are between 10 and 150ms
    durs = offs-ons;
    kills = find(durs>0.15);
    ons(kills)=[];
    offs(kills)=[];
    durs = offs-ons;
    kills = find(durs<0.01);
    ons(kills)=[];
    offs(kills)=[];
    
    keepit=0;
    if (length(ons) > numwind)
        for ii = 1:length(ons)
            p = find(abs(ons(ii:length(ons))-ons(ii))<=wind);
            if (length(p)>=numwind)
                keepit=keepit+1;
            end
        end
        if (keepit>=numnote)
            fprintf(fkeep,'%s\n',fn);
            %disp('keeping...');
        else
            fprintf(fdcrd,'%s\n',fn);
            %disp('discarding...');
        end
    else
        fprintf(fdcrd,'%s\n',fn);
        %disp('discarding...');
    end
end
fclose(fid);fclose(fkeep);fclose(fdcrd);
disp(['done.']);
return
