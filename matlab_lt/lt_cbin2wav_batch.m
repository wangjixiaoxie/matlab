%% lt 10/28/15 - convert batch of song from cbin to wav - puts them in same folder
function lt_cbin2wav_batch(batch, MoveToBirdFolder)


% MoveToBirdFolder=1; then moved up one dir, to a common folder for this
% bird

if ~exist('MoveToBirdFolder', 'var');
    MoveToBirdFolder=0;
end

% =================

fid=fopen(batch);

line=fgetl(fid);

while ischar(line)
    
    if strcmp(line(end-3:end), 'cbin');
       
        if MoveToBirdFolder==0
        cbin2wav(line);
        else
            currdir=pwd;
            slashes=findstr(currdir, '/');
            SaveDir=currdir(1:slashes(end)-1);
            SaveDir2=currdir(slashes(end)+1:end);
            SaveDir=[SaveDir '/SONGS_wav' '/' SaveDir2];
            try
                cd(SaveDir)
                cd(currdir)                
            catch err
                mkdir(SaveDir)
            end
            
            cbin2wav(line, [SaveDir '/' line]);
            disp(['Saved ' line ' to ' SaveDir]);
        end
    end
    
    line = fgetl(fid);
end


