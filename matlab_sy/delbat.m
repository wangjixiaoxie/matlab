% delete files named in batchfile "delfile".  
% written by sam sober, 12/20/05
function [] = delwav(delfile);

if exist(delfile)
    del_fid=fopen(delfile,'r');
end

while 1
    %get filename; check filetype
    soundfile = fscanf(del_fid,'%s',1);
    if isempty(soundfile)
        break
    else
        if isunix
            eval(['!rm ' soundfile]);
        else
            eval(['!del ' soundfile]);
        end
    end
end
fclose(del_fid);
