function lt_make_batch(input,fraction, filetype)
% If input:
% 1: autosort songs based on ampl crossings
% 2: subsample randomly, then autosort (input should be 2,0.2, for 20%)
% 3: write all songs to batch
% 4: make batch with FB hits - i.e. likely real songs
% 5: make batch of all .wav files

if ~exist('filetype', 'var')
    filetype = 'cbin';
end


if ~strcmp(filetype, 'cbin')
    batchname = ['batch_' filetype];
    eval(['!ls *.' filetype ' >  ' batchname]);
else
%makes a batch file
batchname = 'batch';
% batchname = input('What is the name of the batch file?  ', 's');
db_write_batch(batchname)
end

if input==1;
lt_cleandirAuto(batchname,1000,4,4, filetype)
end

if input==2;
    randsamp(batchname,fraction);
    lt_cleandirAuto([batchname '.rand'],1000,4,4);
end

if input==3;
    disp('Batch of all songs made');
end

if input==4;
    lt_rec_files_find_FB_v3(batchname);
end

if input==5;
    wavfilenames=dir('*.wav');
    fid=fopen([batchname '_wav'],'w');
    
    for i=1:length(wavfilenames);
    fprintf(fid,'%s\n',wavfilenames(i).name);
    end
end
    

end

