function [pat]=evgetpat(labels,pattern,outfile);
%[pat]=evgetpat(labels,pattern,outfile);
% pattern matching program
% labels  - string vecotr from songanal which has the labels
% pattern - a string with a grep style pattern to match
% pat  - output pattern match structure  

if (~exist('outfile'))
    outfile = 'TEMP_LABEL_FILE.TEMP';
end
if (exist(outfile,'file'))
    warning(['Will overwrite the label file :',outfile]);
end

fout = fopen(outfile,'w');
if (fout == -1)
    error(['Could not write the output label file!']);
end
fprintf(fout,'%s\n',labels);
fclose(fout);

%run perl scrit to search for regular expressions
eval(['!~/matlab/getpat.pl ',pattern,' ',outfile]);

%search results are in the file : SEARCH_RESULTS.TEMP
%this will load the pattern match results into a structure
fid=fopen('SEARCH_RESULTS.TEMP','r');

%load up the first file num
strind=0;
ln=fgetl(fid);    
pos=findstr(ln,'file number');
while (1)
    if (length(pos)==0)
        break;
    end
    strind = strind + 1;
    fnum = str2num(ln(15:end));
    pat(strind).fname = 'Current Labels File';
    matches={};positions = [];mcnt=0;
    while (1)
        ln=fgetl(fid);
        pos = findstr(ln,'pos');
        if (length(pos)>0)
            positions = [positions;str2num(ln(9:end))];

            ln=fgetl(fid);
            pos = findstr(ln,'match');
            mcnt = mcnt+1;
            match{mcnt} = ln(9:end);
        else
            pat(strind).pos = positions;
            pat(strind).mpats = match;
            pos = findstr(ln,'file number');
            break;
        end
    end
end
fclose(fid);
return;
