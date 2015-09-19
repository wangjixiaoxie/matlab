%this is starting as a general purpose inactivanal script.
%go through each directory, pull out the batch file that you're interested
%in...
%get the first filetime and end file time
%then search the vals struct for files in those times.
%take mean and standard deviation of those values, and n, return value to a
%cell array...

%make the vals struct


%THis path information is specific to the bird
%%it is saved in a summary .mat file in the datasum directory.
%%I need to extract data values.

function [avls graphvals]=inactivanal4(avls,graphvals)


[stimes etimes] = batch2time(avls.pvls,avls.cvl);

    avls=chunkdata2(avls,graphvals)
    avls=calcmeanstd(avls)
    avls=calchist(avls)
    avls=calchistcomb2(avls);
    avls=calcrat2(avls);
    avls=calcmax(avls);
    avls=sortmuruns(avls);
    [avls.wn avls.wnrev]=getwntimes(graphvals.wn,graphvals.wnrev);
   if(isfield(graphvals.wn,'tmlist'))
      [avls.wntms avls.wnptvls]=getwnvls(graphvals); 
   end
    strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)
    strcmd=['save -append ' avls.mtflnm '-analdata.mat avls graphvals' ];
    eval(strcmd);

    function [wntms, ptvls]=getwnvls(graphvals)
        tmlist=graphvals.wn.tmlist;
        ptlist=graphvals.wn.ptlist;
        
        for ii=1:length(tmlist)
           for jj=1:length(tmlist{ii})
                wntms{ii}(jj)=datenum(tmlist{ii}{jj},'yyyy-mm-dd HH') 
                ptvls{ii}(jj)=ptlist{ii}(jj);
           end
        end

    
%plotmeanstdev(chkarrayb, stimes, etimes, graphvals,mnvl,stdv)
% 
% if (combdata==0)
%      strcmd=['cd ' sumpath 'datasum']
%      eval(strcmd);
%      load vals.mat
% end
% 



% 
% histbpre=histc(valsbpre(:,2),edges)
% histbmu=histc(valsmub(:,2),edges)
% histbpost=histc(valsbpost(:,2),edges)
% 
% 
% histbpre=histbpre/sum(histbpre)
% histbmu=histbmu/sum(histbmu)
% histbpost=histbpost/sum(histbpost)
% 
% figure
% stairs(edges,histbpre)
% hold on;stairs(edges,histbmu,'r--')
% 
% figure
% stairs(edges,histbpost)
% hold on;stairs(edges,histbmu,'r--')
% 
% 
% indb=find(valsb(:,2)>3400&valsb(:,2)<4200)
% stdb=std(valsb(indb,2))
% indmb=find(valsmb(:,2)>3400&valsmb(:,2)<4200)
% stdmb=std(valsmb(indmb,2))
