function inds=getmatchind(pat,MAKEBATCH);
% inds=getmatchind(pat,MAKEBATCH);
% 
% finds the places where there were matches
%

if (~exist('MAKEBATCH'))
    MAKEBATCH=0;
end
if (MAKEBATCH==1)
    fid = fopen('TEMPBATCH','w');
end
inds=[];
for ii = 1:length(pat)
    if (length(pat(ii).st)>0)
        inds=[inds;ii];
        if (MAKEBATCH==1)
            ppp=findstr(pat(ii).fn,'.not.mat');
            fprintf(fid,'%s\n',pat(ii).fn(1:ppp(1)-1));
        end
    end
end
if (MAKEBATCH==1)
    fclose(fid);
end

return;