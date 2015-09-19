%This dynamic struct will serve for analysis and plots of 
%rates of consolidation.

% [phsumbs] = mkdynstruct(phsumbs, bs,1:11,'phar')
%type could be 'phar'

function [sumbs] = mkdynstruct(sumbs, bs,ind,type)
ps.numvls=30;
ps.type=type;

for ii=1:length(ind)
    crind=ind(ii);
    [avls]=loadsumdata2(bs, crind);
    [unq_days]=get_uniquedays(avls,ps);
    sumbs(crind).alldays=unq_days-unq_days(1);
    sumbs(crind).basday=unq_days(1);
    outvls=getacvls(sumbs(crind),avls,ps);
    sumbs(crind).initptvls=outvls.initptvls;
    sumbs(crind).endptvls=outvls.endptvls;
    sumbs(crind).init_tms=outvls.init_tms;
    sumbs(crind).end_tms=outvls.end_tms;
    sumbs(crind).initrun=outvls.initrun;
    sumbs(crind).endrun=outvls.endrun;
    
    
    outvls=getmuvls(sumbs(crind),avls,ps);
    sumbs(crind).murun=outvls.murun;
    sumbs(crind).muday=outvls.muday;
    sumbs(crind).muvls=outvls.muvls;
    sumbs(crind).mutms=outvls.mutms;
    sumbs(crind).mutmsinit=outvls.mutmsinit;
end

%returns a list of unique days.
function [unq_days] = get_uniquedays(avls,ps)
    
    if ps.type=='stim'
        
        
        
    %type=pharm    
    else
        dys1=floor(avls.tmvc(:,1));
        dys2=floor(avls.tmvc(:,2));
        comb_dys=[makerow(dys1) makerow(dys2)];
        unq_days=unique(comb_dys);
        unq_days=sort(unq_days);

    end
    
 %input: structure.
 %output: [2 column vec: list of acsf runs for each day, 0 if no run], 
 %identify the earliest and latest acsf run.
 %addvls, tmvls.
 function [outvls] = getacvls(crbs,avls,ps)
     %for each day, identify the earliest and latest acsf run
     alldays=crbs.alldays+crbs.basday;
     numvls=ps.numvls;
     
     
     %one is for beginning of day, two is for end of day
     for outvec=1:2
        for ii=1:length(alldays)
            crday=alldays(ii);
            matchind=find(floor(crbs.tmvc(:,outvec))==crday);
            %now I have a bunch of inds that match this day.
            %find those inds which are acsf runs.
               [outmatch]=find(ismember(matchind,avls.acon));
               matchind=matchind(outmatch); 
            %to get first one.
            [out,ind]=sort(crbs.tmvc(matchind,outvec));
            
            if(~isempty(ind))
                if(outvec==1)
                
                    matchrun=matchind(ind(1));
                   
                    crvls=avls.rawvls{crbs.ntind}{matchrun};
                    outvls.initrun(ii)=matchrun;
                    if(~isempty(crvls))
                    
                        
                        if(length(crvls(:,2))<numvls)
                            ln=length(crvls(:,2))
                        else
                            ln=numvls
                        end
                        outvls.initptvls{ii}=crvls(1:ln,2);
%                         outvls.init_tms{ii}=crvls(1:ln,1)-floor(crvls(1,1));
                        outvls.init_tms{ii}=crvls(1:ln,1);
                    else
                        outvls.initptvls{ii}=[];
                        outvls.init_tms{ii}=[];
                    end
                else
                    matchrun=matchind(ind(end));
                    outvls.endrun(ii) =matchrun;
                    crvls=avls.rawvls{crbs.ntind}{matchrun};
                    if(~isempty(crvls))
                        
                        if(length(crvls(:,2))<numvls)
                            ln=length(crvls(:,2))
                        else
                            ln=numvls
                        end
                    
                        outvls.endptvls{ii}=crvls(end:-1:end-ln+1,2);
%                         outvls.end_tms{ii}=crvls(end:-1:end-ln+1,1)-floor(crvls(end,1));
                        outvls.end_tms{ii}=crvls(end:-1:end-ln+1,1);
                    else
                       outvls.endptvls{ii}=[];
                       outvls.end_tms{ii}=[];
                    end
                    
                    end
            else
                outvls.initptvls{ii}=[];
                outvls.init_tms{ii}=[];
                outvls.endptvls{ii}=[];
                outvls.end_tms{ii}=[];


            end
        end
     end
    
 %Separate function
 %identify stim run(s), muruns on that day, output as an indexed vector.   
 %
 %output mn time, range time of adjvls.    
 function [outvls]=getmuvls(crbs,avls,ps)
     alldays=crbs.alldays+crbs.basday;
%      outvls.murun=[]
%      outvls.muindvl=[]
     ct=1;
     for ii=1:length(alldays)
            crday=alldays(ii);
            matchind=find(floor(crbs.tmvc(:,1))==crday);
            %now I have a bunch of inds that match this day.
            %find those inds which are acsf runs.
               [outmatch]=find(ismember(matchind,avls.muon));
               if(~isempty(matchind(outmatch)))
                    for jj=1:length(outmatch)
                        crvls=avls.adjvls{crbs.ntind}{matchind(outmatch(jj))}
                        if(~isempty(crvls(:,1)))
                            outvls.muvls{ct}=crvls(:,2);
%                             outvls.mutms{ct}=crvls(:,1)-floor(crvls(1,1));
                                outvls.mutms{ct}=crvls(:,1);
                            outvls.mutmsinit(ct)=crvls(1,1)-floor(crvls(1,1));
                            outvls.murun(ct)=[matchind(outmatch(jj))];
                            outvls.muday(ct)=[ crbs.alldays(ii) ]
                            ct=ct+1;
                        
                        end
                    end
               else
%                    outvls.muvls{ct}=[];
%                    outvls.mutms{ct}=[];
%                    outvls.murun(ct)=0;
%                    outvls.muindvl(ct)=0;
%                    ct=ct+1;
                end
               
    end
 
 
     