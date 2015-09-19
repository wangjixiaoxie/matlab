function [asymprunout,cnsec_asympruns,asympvl]=calcrevasympind(aczvls, muzvls,ind,bs,flrwnon,flrwnoff,flrtmvec)

maxdiff=0.75;
%Keep it simple:

%first find all inds that are after the last white noise change.
%then find the median of those inds.
%then find ALL the inds which are within 0.75 sd of that median.

%first, what is the last white noise change.

%indin are the shift runs.
if(~isempty(ind))
tmpind=find(makerow((flrtmvec(ind)>=flrwnon(end)))&makerow(abs(aczvls(ind))<1.5))
possible_asympind=ind(tmpind);

%FOR REV, restrict to thosevalues, less than one from the mean.

%calc median
asympvl=median(aczvls(possible_asympind));
%(this is the asympvl)

%find values where difference is no greater than
%0.75 sd from this median - these are asymp values.

asympdiff=aczvls(possible_asympind)-asympvl;
asympindtmp=find(abs(asympdiff)<0.75);
asympruns=possible_asympind(asympindtmp);
cnsec_asympruns=[];
asymprunout=asympruns;
%provisional mindate
mindate=min(flrtmvec(asympruns));

%are there other runs on this mindate??
%that weren't included.
alldates=flrtmvec(possible_asympind);

% if(~isempty(mindate))
%     tmpind=find(alldates==mindate);
%     samedateind=possible_asympind(tmpind);
% 
% %     if(find(abs(aczvls(samedateind)-asympvl)>0.75))
% %         mindate=mindate+1;
% %     end
% 
%     asympindcortmp=find(flrtmvec(asympruns)>=mindate);
%     asymprunout=asympruns(asympindcortmp);
% 
%     cnsec_asympindtmp=find(flrtmvec(asympruns)>mindate);
%     cnsec_asympruns=asympruns(cnsec_asympindtmp);
% else
%    asympvl=[];
%     asymprunout=[];
%     cnsec_asympruns=[]; 
% end
else
    asympvl=[];
    asymprunout=[];
    cnsec_asympruns=[];
end