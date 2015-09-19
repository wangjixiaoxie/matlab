function [pat]=evgetpatbatch(batchfile,pattern,NOCASE,ADDNOTMAT);
%[pat]=evgetpatbatch(batchfile,pattern,NOCASE,ADDNOTMAT);
% pattern matching program
% batchfile - a string containing the name of a batch file with
%             cbin.not.mat files in it
% pattern   - a string with a grep style pattern to match
% writes a file which with all the labels from all the files in the batch
% file concatenated 

if (~exist('batchfile'))
    batchfile = ['batchnotmat'];
end

if (exist(batchfile,'file'))
    fin = fopen(batchfile,'r');
else
    error(['Batch file does not exist - ',batchfile]);
end
if (fin==-1)
    error(['Could not open batch file - ',batchfile]);
end

if (~exist('ADDNOTMAT'))
    ADDNOTMAT=1;
end
if(~exist('NOCASE'))
    NOCASE=0;
end


Ndatf = 0;
Nmatch=0;
while (1)
    filename = fgetl(fin);
    if (~ischar(filename))
        break;
    end

    if (ADDNOTMAT==1)
        filename=[filename,'.not.mat'];
    end

    if (exist(filename))
        load(filename);
    else
        warning([filename,' does not exist - Skipping to next file!']);
        continue; %skip rest of loop
    end

    if (NOCASE==0)
        [st,en,tmp,matstr]=regexp(labels,pattern);
    else
        [st,en,tmp,matstr]=regexpi(labels,pattern);
    end
    
    Nmatch=Nmatch+length(st);
    
    Ndatf = Ndatf + 1;
    pat(Ndatf).fn     = filename;
    pat(Ndatf).labels = labels;
    pat(Ndatf).match  = matstr;
    pat(Ndatf).start  = st;
    pat(Ndatf).end    = en;
end
disp(['Found ',num2str(Nmatch),' matches']);
fclose(fin);
return;
