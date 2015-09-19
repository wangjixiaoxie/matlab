function [aclist,mulist]=groupadj(timelist,acind,muind)


%groupadj(timelist,acind,muind)
%this code will go through each all the acsf start_times and find the unique
%days to create


%for all the acsf inds
%find the unique days.
actimes=timelist(acind,1)
fl_actimes=floor(actimes);
fl_alltimes=floor(timelist);
acdy=unique(fl_actimes);
acindlst=[];
for ii=1:length(acdy)
    actmplst=find(fl_actimes==acdy(ii))
    
    %sort the times in ascending order;
    [dm, sortind]=sort(actimes(actmplst));
    actmplst=acind(actmplst(sortind));
    if (length(actmplst)>2)
        actmp=[actmplst(1:2);actmplst(2:3)]
    
    elseif length(actmplst)==1
          actmp=[actmplst 0]  
    else
            actmp=actmplst(1:2);
    end
            
     
    
    acindlst=[acindlst;actmp]
end

%these are the acsf inds...now search for mu list, which are greater, and
%which floor is same day (if more than one pick the lower one)...this is
%the pair.
muindlst=[];
mutimes=timelist(muind,1);
fl_mutimes=floor(mutimes);
for ii=1:length(acindlst(:,1))
    indvl=acindlst(ii,1);
    muindvl=find(mutimes>timelist(indvl,1) &(fl_mutimes==fl_alltimes(indvl,1)))
   if isempty(muindvl)
       mutmp=0;
   else
       mutmp=muind(muindvl);
   end
    muindlst=[muindlst mutmp(1)];
end
aclist=acindlst;
mulist=muindlst;

