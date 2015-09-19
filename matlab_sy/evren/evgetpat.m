function [pat]=evgetpat(labelstruct,pattern,NOCASE);
%[pat]=evgetpat(labelstruct,pattern,NOCASE);
% pattern matching program
% labels  - string vecotr from songanal which has the labels
% pattern - a string with a grep style pattern to match
% pat  - output pattern match structure  

if (~exist('NOCASE'))
    NOCASE=0;
end
Nmatch=0;
for ii = 1:length(labelstruct)
    lbl=labelstruct(ii).labels;
    if (NOCASE==0)
        [st,en,tmp,mat]=regexp(lbl,pattern);
    else
        [st,en,tmp,mat]=regexpi(lbl,pattern);
    end
    pat(ii).fn=labelstruct(ii).fn;
    pat(ii).st=st;
    pat(ii).en=en;
    pat(ii).match=mat;
    Nmatch=Nmatch+length(st);
end
disp(['n match = ',num2str(Nmatch)]);
return;
