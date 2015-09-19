%CALL THIS FUNCTION WITH SPECIFIC INPUTS set in a bird-specific .m file
%Data saved to a .mat file
function []=synanal3(avls)
%THIS IS A HACK, specifying which note to analyze for pitch analysis
% PTNT=3;


MKFV=avls.MKFV;
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
        
        if(avls.SEC(ii)==0)
            vls{crnt}=getvals(fvpt{crnt},1,'TRIG')
        else
            vls{crnt}=getvals_sec(fvpt{crnt},1,'TRIG')
        end
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
combtms=[]
combvls=[];
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
    ptvls=[vls{1} ;vls{avls.ANALPTNT}]
    [out,ind]=sort(ptvls(:,1))
    combvls=[combvls; ptvls(ind,:)];
end


%now set up the individual responsevecs for syntax calculations.
%divide into days
unq_days=unique(floor(combtms));
for ii=1:length(unq_days)
   res_vec=[];
    crday=unq_days(ii);
   crdayind=find(floor(combtms)==crday)
   ptvlsind=find(floor(combvls(:,1))==crday)
   res_vec=zeros(length(crdayind),1);
   for ntind=1:length(avls.NT)
       
        indnt{ntind}=find(comblabels(crdayind)==avls.NT{ntind});
        outnotect(ii,ntind)=length(indnt{ntind})
   end
  
   res_vec(indnt{avls.SEQNT})=1;
   dayout(ii)=crday;
   if(avls.DOCONT_ANAL)
        runanalysis(res_vec,1,0.5);
        cmd=['load resultsindividual.mat'];
        eval(cmd);
        vlsout{ii}=combvls(ptvlsind,:);
  
        tmscomb=[combtms(crdayind);combtms(crdayind(end))];
   %average values with the same "time", i.e. in same song    
        unqtms=unique(tmscomb);
        for tmind=1:length(unqtms)
            selind=find(tmscomb==unqtms(tmind));
            tms{ii}(tmind)=unqtms(tmind);
            pout{ii}(tmind)=mean(p(selind));
            p05out{ii}(tmind)=mean(p05(selind));
            p95out{ii}(tmind)=mean(p95(selind));
        end
   end
end

%provide indices for bas and rec ind
for ii=1:length(avls.bastms)
   basbnds=datenum(avls.bastms{ii});
   wnbnds=datenum(avls.wntms{ii});
   recbnds=datenum(avls.rectms{ii});
   
    basind{ii}=find(dayout>=basbnds(1)&dayout<=basbnds(2));
   wnind{ii}=find(dayout>=wnbnds(1)&dayout<=wnbnds(2));
   recind{ii}=find(dayout>=recbnds(1)&dayout<=recbnds(2));
end



cmd=['cd ' avls.baspath ]
eval(cmd);
if(avls.DOCONT_ANAL)
    cmd=['save -append  sumdata.mat  tms pout p05out p95out vlsout outnotect dayout'];

    eval(cmd);
else
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat  outnotect dayout'];
        eval(cmd); 
    else
        cmd=['save  sumdata.mat  outnotect dayout'];
        eval(cmd); 
    end
end
    







