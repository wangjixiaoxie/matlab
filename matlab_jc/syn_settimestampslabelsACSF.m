%CALL THIS FUNCTION WITH SPECIFIC INPUTS set in a bird-specific .m file
%Data saved to a .mat file
function []=syn_settimestampslabelsACSF(avls)
MKFV=avls.MKSYN;
for ii=1:length(MKFV)
    crind=MKFV(ii);
    crpt=avls.pvls{crind}
    crbt=avls.cvl{crind}
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath crpt]
    else
        cmd=['cd ' crpt]
    end
    eval (cmd);      labels=notestats_tw(crbt);
    %initialize timestamps to zero
    alltms=zeros(length(labels),1);
    clear outstruct 
    for crnt=1:length(avls.NT)
        NT=avls.NT{crnt};PRENT='';PSTNT='';
        crtmshft=avls.tbinshft(crnt)
        crbt=avls.cvl{crind};
        crfbins=avls.fbins{crnt};
        FFT=avls.NFFT(crnt);
        %structure for pitch measurments
        fvpt{crnt}=findwnote4(crbt,NT, PRENT, PSTNT, crtmshft,crfbins, FFT,1,'obs0');        
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

% For each folder, calculate the learning curve
for ii=1:length(avls.pvls)
    res_vec=[];
    % load data
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
    % outstruct contains labels and times
        comblabels=outstruct{1};
        tmscomb=outstruct{2};
   res_vec=zeros(length(tmscomb),1);
   % for each targeted note == TRUE
       for i=1:length(avls.targnts)
        indnt=find(comblabels==avls.targnts(i));
        res_vec(indnt)=1;
       end
   if ii==1
       pnull=0.5;
   else
       pnull=mean(p(round(end/2):end)); % could be improved
   end
   runanalysis(res_vec,1,pnull,0.005);
   cmd=['load resultsindividual.mat'];
   eval(cmd);
   %average values with the same "time", i.e. in same song    
   tms{ii}=tmscomb;
   pout{ii}=p(1:end-1);
   p05out{ii}=p05(1:end-1);
   p95out{ii}=p95(1:end-1);
%    unqtms=unique(tmscomb);
%    for tmind=1:length(unqtms)
%        selind=find(tmscomb==unqtms(tmind));
%        tms{ii}(tmind)=unqtms(tmind);
%        pout{ii}(tmind)=mean(p(selind));
%        p05out{ii}(tmind)=mean(p05(selind));
%        p95out{ii}(tmind)=mean(p95(selind));
%    end
end

cmd=['mkdir ' avls.baspath 'datsum']
eval(cmd);
cmd=['cd ' avls.baspath 'datsum']
eval(cmd);
cmd=['save  sumdataACSF.mat  tms pout p05out p95out '];
eval(cmd);


   
    







