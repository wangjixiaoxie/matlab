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

function [avls graphvals]=inactivanal(avls,graphvals,fvcomb);

strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)

%this gives you access to fvcombs.
strcmd=['load ' avls.mtflnm '-analdata.mat '];
    eval(strcmd);
    

[stimes etimes] = batch2time(avls.pvls,avls.cvl);

%modified 10/16/07, load each batch file now and find the meanval for that
%batch file.


    
    ofst=avls.offset/24;
    acofst=avls.acoffset/24
    
    %fvcomb file already loaded
    %if there is a timon for the file, throw out data until the offset
    for ii=1:length(avls.analind)
        btind=avls.analind(ii)
        clear fv
        strcmd=['cd '  avls.pvls{btind}];
        eval(strcmd);
        strcmd=['load ' avls.pvls{btind} avls.cvl{btind} '.mat fv']
        eval(strcmd);
        for jj=1:length(avls.NT)
            vals=getvals(fv{jj},1,'TRIG');
            if(graphvals.timon{ii})
                t_strt=datenum([graphvals.date{ii} ' ' graphvals.timon{ii}],'yyyy-mm-dd  HH:MM:SS')
                t_end=datenum([graphvals.date{ii} ' ' graphvals.timoff{ii}],'yyyy-mm-dd  HH:MM:SS')
                vlsind=find(vals(:,1)>(t_strt+ofst));
                 chkarrayadj{jj}{btind}=vals(vlsind,:);
            elseif((graphvals.acon(ii)))
                acind=graphvals.acon(ii);
                t_strt=datenum([graphvals.date{acind} ' ' graphvals.timoff{acind}], 'yyyy-mm-dd HH:MM:SS')
                vlsind=find(vals(:,1)>(t_strt+acofst));
                chkarrayadj{jj}{btind}=vals(vlsind,:);
            
            else        
                    
               chkarrayadj{jj}{btind}=vals;
            end
        end
    end
        
    
    
    
    for ii=1:length(avls.analind)
        btind=avls.analind(ii);
        clear fv
        strcmd=['cd '  avls.pvls{btind}];
        eval(strcmd);
        strcmd=['load ' avls.pvls{btind} avls.cvl{btind} '.mat fv']
        eval(strcmd);
        for jj=1:length(avls.NT)
            chkarray{jj}{btind}=getvals(fv{jj},1,'TRIG')
        end
    end
    
    
    


    %first note and second note, need to generalize

    %get mean and standard deviation of each time
    clear mnvl stdv
for noteind=1:length(chkarray)
    for jj=1:length(avls.analind) 
        ssnind=avls.analind(jj)
        
        %clear out the outliers
        vals=chkarrayadj{noteind}{ssnind}(:,2);
        valsadjind=find(vals>avls.fbins{noteind}(1)&vals<avls.fbins{noteind}(2));
        valsmod=vals(valsadjind);
        mnvl{noteind}(ssnind)=mean(valsmod)
        stdv{noteind}(ssnind)=std(valsmod)
    
        hstout{noteind}{ssnind}=histc(valsmod,graphvals.edges{noteind})
        nnote{noteind}{ssnind}=length(valsmod)
    
    end
end


%make histogram, using graphvals.edges



clear graphvals.colind
graphvals.colind=graphvals.colvals;
for ii=1:graphvals.numcol
    for jj=1:length(graphvals.colind{ii})
        ind=graphvals.colind{ii}(jj)
        graphvals.collist(ind)=graphvals.col(ii);   
      end
end
avls.chkarray=chkarray
avls.mnvl=mnvl
avls.hstout=hstout
avls.stdv=stdv
avls.stimes=stimes
avls.etimes=etimes
avls.nnote=nnote

% if(graphvals.chunkdata)
%     avls.chunkout=chunkdata(avls,graphvals)
% end



strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)
    strcmd=['save -append ' avls.mtflnm '-analdata.mat avls graphvals' ];
    eval(strcmd);

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
