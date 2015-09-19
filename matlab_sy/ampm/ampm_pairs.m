%piloting this just for bk61w42...  

%mornevescatter

%find the matching am/pm pairs
%do this by going for each am and checking for the appropriate pm.

function [am_pmout,am_muout]=ampm_pairs (phsumbs,bs,bsind,ps)

%foram_pmanal.
 [am_pmout]=findampm_pairs(phsumbs,bsind);
 
 %foram_muanal
 [am_muout]=findam_mupairs(phsumbs,bsind);


 
function [am_pmout]=findampm_pairs(phsumbs,bsind)
    for ii=1:length(bsind)
        crbsind=bsind(ii);
        crbs=phsumbs(crbsind);
        
        %outruns{1} is shiftruns
        %outruns{2} is asympruns
        [outruns]=getruns(crbs)
        for ii=1:2
            cr_runs=outruns{ii};
            [amruns,pmruns]=getampmruns(cr_runs,crbs);
            %asymp
            if(ii==1)
                am_pmout.shift.amruns=amruns;
                am_pmout.shift.pmruns=pmruns;
            %shift
            else
                am_pmout.asymp.amruns=amruns;
                am_pmout.asymp.pmruns=pmruns;
            end
        end
    end

function [am_muout]=findam_mupairs(phsumbs,bsind)
    for ii=1:length(bsind)
        crbsind=bsind(ii);
        crbs=phsumbs(crbsind);
        
        %outruns{1} is shiftruns
        %outruns{2} is asympruns
        [outruns]=getruns(crbs)
        for ii=1:2
            cr_runs=outruns{ii};
            [amruns,muruns]=getammuvls(cr_runs,crbs);
            %asymp
            if(ii==1)
            
            am_muout.shift.amruns=amruns;
            am_muout.shift.muruns=muruns;
            %shift
            else
                am_muout.asymp.amruns=amruns;
                am_muout.asymp.muruns=muruns; 
                
            end
        end   
    end


function [outruns]=getruns(crbs)
%get only shift runs including first asymptote run.
outruns{1}=[];
outruns{2}=[];

for ii=1:length(crbs.allasympruns)
    crasymp=crbs.allasympruns{ii};
    outruns{2}=[outruns{2} ;crasymp]
    
    
    mulistind=find(ismember(crbs.mulist,crasymp));
    if(~isempty(mulistind))
        mulim(ii)=mulistind(1);
    else
        mulim(ii)=0;
    end
        
end
for ii=1:length(crbs.allshiftruns)
    crshift=crbs.allshiftruns{ii};
    
    %limit to first asymptote run
    if((mulim(ii)))
%        acinitruns=crbs.aclist(1:mulim(ii),1);
%        acfinruns=crbs.aclist(1:mulim(ii),2);
%        muruns=crbs.mulist(1:mulim(ii));
%        nozeroac=find(acfinruns~=0);
%        nozeromu=find(muruns~=0);
%        acfinruns=acfinruns(nozeroac);
%        mufinruns=muruns(nozeromu);
         tmvc=crbs.flrtmvec;
         tmbnd=tmvc(crbs.mulist(mulim(ii)));
         
        crshiftind=find(tmvc(crshift)<=tmbnd);
       outruns{1}=[outruns{1} crshift(crshiftind)]
    end
end
        
function [amruns,pmruns]=getampmruns(cr_runs,crbs)
%first find the am values.
amrunsind=find(ismember(crbs.aclist(:,1),cr_runs));
amruns=crbs.aclist(amrunsind,1);
if(isempty(amruns))
    pmruns=0
end

%now look for a pm pair from day before.
%for each of the amruns, find the aclistind
for ii=1:length(amruns)
   aclistind=find(crbs.aclist(:,1)==amruns(ii))
   if(aclistind-1)
        pmlistind=aclistind-1;
        pmruns(ii)=crbs.aclist(pmlistind,2);
        if(~(pmruns(ii)))
            pmruns(ii)=crbs.aclist(pmlistind,1)
        end
   else
       pmruns(ii)=0;


   end
end

function [amruns,muruns]=getammuvls(cr_runs,crbs)
amrunsind=find(ismember(crbs.aclist(:,1),cr_runs));
amruns=crbs.aclist(amrunsind,1);
if(isempty(amruns))
    muruns=[];
end
for ii=1:length(amruns)
    murunsvl=crbs.mulist(amrunsind(ii));
    if(~isempty(murunsvl))
        muruns(ii)=murunsvl;
    else
        muruns(ii)=0
        
    end
end




    %then plot the values with standard error


%for am mu scatter
%loop through the mu and find the matching am.