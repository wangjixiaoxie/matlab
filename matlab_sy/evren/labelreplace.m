function [labelstructout]=labelreplace(labelstruct,pattern,newpattern,DOLABEL,NOCASE);
%[labelstructout]=labelreplace(labelstruct,pattern,newpattern,DOLABEL,NOCASE);

if (~exist('NOCASE'))
    NOCASE=0;
end

if (~exist('DOLABEL'))
    DOLABEL=0;
end

for ii = 1:length(labelstruct)
    lbl=labelstruct(ii).labels;
    
    if (NOCASE==0)
        [st,en,tmp,mat]=regexp(lbl,pattern);
    else
        [st,en,tmp,mat]=regexpi(lbl,pattern);
    end
    
    for jj=1:length(st)
        if (length(newpattern)==(en(jj)-st(jj)+1))
            lbl(st(jj):en(jj))=newpattern;
        else
            disp(['Wrong size? fn : ',labelstruct(ii).fn,' mat : ',mat{jj},...
                                                      ' newpat : ',newpattern]);
        end
    end
    
    labelstruct(ii).labels=lbl;
    labels=lbl;
    if (DOLABEL==1)
        if (exist(labelstruct(ii).fn,'file'))
            cmd=['save -append ',labelstruct(ii).fn,' labels'];
            eval(cmd);
        end
    end
end

labelstructout=labelstruct;
return; 