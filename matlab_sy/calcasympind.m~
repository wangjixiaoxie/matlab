function [asympindout,cnsec_asympind,asympvl]=calcasympind(avls,aczvls,ps )

ntind=ps.ntind;
indin=ps.indin;
flrwnon=ps.flrwnon;
flrwnoff=ps.flrwnoff;
flrtmvec=ps.flrtmvec;

TYPE=ps.type;


%NEED TO FIX THIS by creating a list of possibleacind only, but then adding
%automatically any mu runs in between...
if(TYPE=='phar')
   indin=find(ismember(avls.aclist,indin))
   indin=avls.aclist(indin);
   aczvls=avls.allz(:,ps.ntind);

end
maxdiff=1;
%Keep it simple:

%first find all inds that are after the last white noise change.
%then find the median of those inds.
%then find ALL the inds which are within 0.75 sd of that median.

%first, what is the last white noise change.

%indin are the shift runs.

tmpind=find(flrtmvec(indin)>=flrwnon(end))
possible_asympind=indin(tmpind);

%calc median
asympvl=nanmedian(aczvls(possible_asympind));
%(this is the asympvl)

%find values where difference is no greater than
%0.75 sd from this median - these are asymp values.

asympdiff=aczvls(possible_asympind)-asympvl;
asympindtmp=find(abs(asympdiff)<maxdiff);
asympind=possible_asympind(asympindtmp);

%provisional mindate
mindate=min(flrtmvec(asympind));

%are there other runs on this mindate??
%that weren't included.
%RESTRICT TO FIRST DATE WHEN ALL RUNS WERE WITHIN BOUNDS
alldates=flrtmvec(possible_asympind);
if(~isempty(mindate))
    tmpind=find(alldates==mindate);
    samedateind=possible_asympind(tmpind);
    if(find(abs(aczvls(samedateind)-asympvl)>maxdiff))
        mindate=mindate+1;
    end

    asympindcortmp=find(flrtmvec(asympind)>=mindate);
    asympindout=asympind(asympindcortmp);
    if(TYPE=='phar')
       munozero=find(avls.mulist>0);
       muvls=avls.mulist(munozero);
       minactm=min(floor(avls.tmvc(asympindout,1)));
       maxactm=max(floor(avls.tmvc(asympindout,1)));
       muind=find((floor((avls.tmvc(muvls,1)))>=minactm)&(floor(avls.tmvc(muvls,1))<=maxactm))
       asympindout=[makerow(asympindout) makerow(muvls(muind))];
    end

    cnsec_asympindtmp=find(flrtmvec(asympindout)>mindate);
    cnsec_asympind=asympindout(cnsec_asympindtmp);
else
    asympindout=[];
    cnsec_asympind=[];
    asympvl=[];
end
tst=1;