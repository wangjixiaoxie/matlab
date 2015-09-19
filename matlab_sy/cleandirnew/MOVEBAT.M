% moves files named in batchfile "delfile" to directory in move_path  
% written by sam sober, 12/20/05
function [] = delwav(delfile,move_path);

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
error
            %            eval(['!rm ' soundfile]);
        else
            if isdir(move_path)
            eval(['!move ' soundfile ' ' move_path]);
            else
                error([move_path ' is not a valid path'])
            end
        end
    end
end
fclose(del_fid);
