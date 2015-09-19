function [pat]=evgetpatbatch(batchfile,pattern,outfile);
%[pat]=evgetpatbatch(batchfile,pattern,outfile);
% pattern matching program
% batchfile - a string containing the name of a batch file with
%             cbin.not.mat files in it
% pattern   - a string with a grep style pattern to match
% writes a file which with all the labels from all the files in the batch
% file concatenated 

if (~exist('batchfile'))
    batchfile = ['batchnotmat'];
end

if (~exist('outfile'))
    outfile = 'TEMP_LABEL_FILE.TEMP';
end
if (exist(outfile,'file'))
    warning(['Will overwrite the label file :',outfile]);
end

if (exist(batchfile,'file'))
    fin = fopen(batchfile,'r');
else
    error(['Batch file does not exist - ',batchfile]);
end
if (fin==-1)
    error(['Could not open batch file - ',batchfile]);
end

fout = fopen(outfile,'w');
if (fout == -1)
    error(['Could not write the output label file!']);
end

Ndatf = 0;
while (1)
    filename = fgetl(fin);
    if (~ischar(filename))
        break;
    end
    
    if (exist(filename))
        load(filename);
    else
        warning([filename,' does not exist - Skipping to next file!']);
        continue; %skip rest of loop
    end
    
    fprintf(fout,'%s\n',labels);
    Ndatf = Ndatf + 1;
    datafiles{Ndatf,1} = filename; %cell arrray of data file names
end
fclose(fin);fclose(fout);

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
    pat(strind).fname = datafiles{fnum};
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
