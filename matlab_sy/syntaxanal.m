function [synstruct] = syntaxanal(fv,translist)
    daylist=[];
   for ii=1:length(fv)
        dvl=fn2datenum(fv(ii).fn);
        fldvl=floor(dvl);
        daylist=[daylist fldvl];
    end
    vlot=unique(daylist);
    

    %loop through all the days
for ii=1:length(vlot)
    
    ind=find(daylist==(vlot(ii)));
    %loop through all the indices for that day
    for kk=1:length(translist)
            cmd=['synstruct{ii}.' translist{kk} '=0'];
            eval(cmd);
    end
    
    for jj=1:length(ind)
        indvl=ind(jj);
        %for each index check the transition.
        
        for kk=1:length(translist)
            if(fv(indvl).trans(1)==translist{kk})
                cmd=['synstruct{ii}.' translist{kk} '=synstruct{ii}.' translist{kk} '+1'];
                eval(cmd)
            end
        end
    synstruct{ii}.day=vlot(ii)
    synstruct{ii}.translist=translist;
    synstruct{ii}.ntrans=length(ind);
    end
end
