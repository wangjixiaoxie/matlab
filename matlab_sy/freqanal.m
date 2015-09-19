function [mnvl,stdv,htrt,vlot]= freqanal(fv)
    daylist=[];
   for ii=1:length(fv)
        dvl=fn2datenum(fv(ii).fn);
        fldvl=floor(dvl);
        daylist=[daylist fldvl];
    end
    vlot=unique(daylist);
    vals=getvals(fv,1,'TRIG');
    for ii=1:length(vlot)
    
        ind=find(daylist==(vlot(ii)));
    
        
            
            mnvl(ii)=mean(vals(ind,2))
            stdv(ii)=std(vals(ind,2))
            nnotes(ii)=length(ind)
            if(isempty(find(vals(ind,3)==1)))
                indht=0;
            else
                indht=find(vals(ind,3)==1)
            end
                htrt(ii)=length(indht)/nnotes(ii);
        
            
    end