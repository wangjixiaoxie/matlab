%CALL THIS FUNCTION WITH SPECIFIC INPUTS set in a bird-specific .m file
%Data saved to a .mat file
function []=syn_settimestampslabelsCURVE(avls)

%%%% Part 1 - extract data
            MKFV=avls.MKSYN;
            for ii=1:length(MKFV)
                crind=MKFV(ii);
                crpt=avls.pvls{crind} % path
                crbt=avls.cvl{crind} % batch
                if(isfield(avls,'baspath'))
                    cmd=['cd ' avls.baspath crpt]
                else
                    cmd=['cd ' crpt]
                end
                eval (cmd);     
                labels=notestats_tw(crbt);
                %initialize timestamps to zero
                alltms=zeros(length(labels),1);
                clear outstruct 
                for crnt=1:length(avls.NT)
                    NT=avls.NT{crnt};PRENT='';PSTNT='';
                    crbt=avls.cvl{crind};
                    %structure for pitch measurments
                    fvpt{crnt}=findwnote4(crbt,NT, PRENT, PSTNT,0,[2000 3000],256,1,'obs0');        
                    %need to get a time stamp for all these different n     
                    vls{crnt}=getvals(fvpt{crnt},1,'TRIG')
                    %indices on labels which are this note.
                    labelind=find(ismember(labels,NT));
                    alltms(labelind)=vls{crnt}(:,1);
                end
                %save these values to a file
                %save outstruct,fvpt,vls
                indnonzero=find(alltms>0);
                outstruct{1}=labels(indnonzero);
                outstruct{2}=alltms(indnonzero);

                  cmd=['save  ' crbt '.mat outstruct fvpt vls'];
                  eval(cmd);
            end
%combine data into one master struct.
%set up a series of pvalues to plot.
comblabels=[];
combtms=[];
for ii=1:length(avls.pvls)
    
    %load the batch file.
    crpt=avls.pvls{ii}
    crbt=avls.cvl{ii}
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath crpt]
    else
        cmd=['cd ' crpt]
    end
    eval (cmd);  
    
    cmd =['load ' crbt '.mat']
    eval(cmd);
    comblabels=[comblabels;outstruct{1}]
    combtms=[combtms;outstruct{2}]
end

%%%% Part 2: For each day , calculate the learning curve
unq_days=unique(floor(combtms));
for ii=1:length(unq_days)
   res_vec=[];
    crday=unq_days(ii);
   crdayind=find(floor(combtms)==crday)
   res_vec=zeros(length(crdayind),1)
   dayout{ii}=datevec(crday)
   % for each targeted note == TRUE
   for i=1:length(avls.targnts)
       indnt=find(comblabels==avls.targnts(i));
       res_vec(indnt)=1;
   end
   if ii==1
       pnull=0.5;
   else
       pnull=mean(p(round(end/4):end)); % could be improved
   end
   runanalysis(res_vec,1,pnull,0.005);
   cmd=['load resultsindividual.mat'];
   eval(cmd);
     tmscomb=[combtms(crdayind);combtms(crdayind(end))]; 
%    pout{ii}=p;
%    p05out{ii}=p05;
%    p95out{ii}=p95;
%    
% %         % no averaging
%             tms{ii}=tmscomb;
%             pout{ii}=p(1:end-1);
%             p05out{ii}=p05(1:end-1);
%             p95out{ii}=p95(1:end-1);
%             pout=pout';
%             p05out=p05out';
%             p95out=p95out';
   tmscomb=[combtms(crdayind);combtms(crdayind(end))];
%    average values with the same "time", i.e. in same song    
   unqtms=unique(tmscomb);
   for tmind=1:length(unqtms)
       selind=find(tmscomb==unqtms(tmind));
       tms{ii}(tmind)=unqtms(tmind);
       pout{ii}(tmind)=mean(p(selind));
       p05out{ii}(tmind)=mean(p05(selind));
       p95out{ii}(tmind)=mean(p95(selind));
   end
end
% 
%         for ii=1:length(avls.pvls)
%             res_vec=[];
%             % load data
%             crpt=avls.pvls{ii}
%             crbt=avls.cvl{ii}
%             if(isfield(avls,'baspath'))
%                 cmd=['cd ' avls.baspath crpt]
%             else
%                 cmd=['cd ' crpt]
%             end
%             eval (cmd);
% 
%             cmd =['load ' crbt '.mat']
%             eval(cmd);
%             % outstruct contains labels and times
%             comblabels=outstruct{1};
%             tmscomb=outstruct{2};
%             res_vec=zeros(length(tmscomb),1);
%             % for each targeted note == TRUE
%             for i=1:length(avls.targnts)
%                 indnt=find(comblabels==avls.targnts(i));
%                 res_vec(indnt)=1;
%             end
%             if ii==1
%                 pnull=0.5;
%             else
%                 pnull=mean(p(round(end/2):end)); % could be improved
%             end
%             runanalysis(res_vec,1,pnull,0.005);
%             cmd=['load resultsindividual.mat'];
%             eval(cmd);
%             
%         % no averaging
%             tms{ii}=tmscomb;
%             pout{ii}=p(1:end-1);
%             p05out{ii}=p05(1:end-1);
%             p95out{ii}=p95(1:end-1);
%             pout=pout';
%             p05out=p05out';
%             p95out=p95out';
%             
%         %average values with the same "time", i.e. in same song
%             %    unqtms=unique(tmscomb);
%             %    for tmind=1:length(unqtms)
%             %        selind=find(tmscomb==unqtms(tmind));
%             %        tms{ii}(tmind)=unqtms(tmind);
%             %        pout{ii}(tmind)=mean(p(selind));
%             %        p05out{ii}(tmind)=mean(p05(selind));
%             %        p95out{ii}(tmind)=mean(p95(selind));
%             %    end
%         end
        
%%%% Part 3: Save the data
        cmd=['mkdir ' avls.baspath 'datsum']
        eval(cmd);
        cmd=['cd ' avls.baspath 'datsum']
        eval(cmd);
        cmd=['save  sumdataCURVES.mat  tms pout p05out p95out '];
        eval(cmd);


   
    







