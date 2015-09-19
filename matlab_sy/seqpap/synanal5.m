%CALL THIS FUNCTION WITH SPECIFIC INPUTS set in a bird-specific .m file
%Data saved to a .mat file
function []=synanal2(avls)
%THIS IS A HACK, specifying which note to analyze for pitch analysis
% PTNT=3;


MKFV=avls.MKFV;
PTIND=avls.PTIND
PTANAL=intersect(MKFV,PTIND);

%first just do pitch analysis
for ii=1:length(PTANAL)
    %first do just those PTIND
    %write this to it's own data struct
    
    crind=PTANAL(ii);
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
        
    
        NT=avls.NT{crnt};PRENT=avls.PRENT{crnt};PSTNT='';
        crtmshft=avls.tbinshft(crnt)
        crbt=avls.cvl{crind};
        crfbins=avls.fbins{crnt};
        FFT=avls.NFFT(crnt);
        %structure for pitch measurments
        fvpt{crnt}=findwnote4(crbt,NT, PRENT, PSTNT, crtmshft,crfbins, FFT,1,'obs0');
        
        %this controls which version of get vals is used.
        
        if(avls.EV4(ii)==0)
            vls{crnt}=getvals(fvpt{crnt},1,'TRIG')
        else
            vls{crnt}=getvals_sec(fvpt{crnt},1,'TRIG')
        end
        %indices on labels which are this note.
        labelind=find(ismember(labels,NT));
        if(~isempty(PRENT))
            outind=find(labels(labelind-1)==PRENT)
            labelind=labelind(outind);
        
       
        end
        if(~isempty(labelind))
            if(~isempty(vls{crnt}))
                if(length(labelind)==length(vls{crnt}(:,1)))
                alltms(labelind)=vls{crnt}(:,1);
                end
                
            end
         end
        
    end
    %save these values to a file
    %save outstruct,fvpt,vls
    indnonzero=find(alltms>0);
    outstruct{1}=labels(indnonzero);
    outstruct{2}=alltms(indnonzero);
    
      cmd=['save  ' crbt '.mat outstruct fvpt vls'];
      eval(cmd);
end
 

%then do just sequence analysis
SEQIND=avls.SEQIND;
SEQANAL=intersect(SEQIND,MKFV)
    
%combine data into one master struct.
%set up a series of pvalues to plot.
comblabels=[];
combtms=[]
combvls=[];
if(~isempty(SEQANAL))
for ii=1:length(SEQANAL)
    crind=SEQANAL(ii);
    %load the batch file.
    crpt=avls.pvls{crind}
    crbt=avls.cvl{crind}
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
    ptvls=[vls{avls.ANALPTNT} ;vls{avls.ANALPTNT}]
    [out,ind]=sort(ptvls(:,1))
    combvls=[combvls; ptvls(ind,:)];
end


%now set up the individual responsevecs for syntax calculations.
%divide into daysed
if (avls.SPLITDAYS==0)
    unq_days=unique(floor(combtms));
else
    unq_days=unique(floor(combtms/.5));
end
if(~isempty(unq_days))
for ii=1:length(unq_days)
    if(ii==1)
        SEEDMEAN=0.5;
    end
   res_vec=[];
    crday=unq_days(ii);
   if(avls.SPLITDAYS)
       crdayind=find(floor(combtms/.5)==crday)
       ptvlsind=find(floor(combvls(:,1)/.5)==crday)
   else
       crdayind=find(floor(combtms)==crday)
       ptvlsind=find(floor(combvls(:,1))==crday)
   end
    vlsout{ii}=combvls(ptvlsind,:);
   res_vec=zeros(length(crdayind),1);
   for ntind=1:length(avls.NT)
       
            indnt{ntind}=find(comblabels(crdayind)==avls.NT{ntind});
            outnotect(ii,ntind)=length(indnt{ntind})
       
           indnt{ntind}=find(comblabels(crdayind)==avls.NT{ntind});
            outnotect(ii,ntind)=length(indnt{ntind})
   end
  
   res_vec(indnt{avls.SEQTRGNT})=1;
   if(avls.SPLITDAYS==0)
       dayout(ii)=crday;
   else
       dayout(ii)=crday/2;
   end
       
       if(avls.DOCONTSYNANAL||avls.DOCONTPTANAL)
        %rhis is the sequence analysis
        %res_vec is a series of zeros and ones; need to change expected
        %probability from 0.5 to prob. from prior day
        if(avls.DOCONTSYNANAL)
            if(~isempty(find(res_vec>0)))
                runanalysis(res_vec,1,SEEDMEAN);
                cmd=['load resultsindividual.mat'];
                eval(cmd);
            
                tmscomb=[combtms(crdayind);combtms(crdayind(end))];
   %average values with the same "time", i.e. in same song    
                unqtms=unique(tmscomb);
                for tmind=1:length(unqtms)
                    selind=find(tmscomb==unqtms(tmind));
                    cntseq(ii).tms(tmind)=unqtms(tmind);
                    cntseq(ii).svl(tmind)=mean(p(selind));
                    cntseq(ii).s05(tmind)=mean(p05(selind));
                    cntseq(ii).s95(tmind)=mean(p95(selind));
                end
            SEEDMEAN=mean(cntseq(ii).svl);
            else
                SEEDMEAN=0.5;
            end
        end
        %this is the continuous pitch analysis.
        if(avls.DOCONTPTANAL)
            WIN=15;
            if(~isempty(vlsout{ii}))
                if(length(vlsout{ii}(:,1))>3*WIN)
%                     cntpt(ii).tms=runave2(vlsout{ii}(:,1),WIN);
                      startind=ceil(WIN/2);
                      cntpt(ii).tms=vlsout{ii}(startind:(end-ceil(WIN/2)));
%                   cntpt(ii).fmed=runave2(vlsout{ii}(:,2));
                    [datout]=calcerbndscnt(vlsout{ii}(:,2),WIN);
                    if(~isempty(datout))
                        cntpt(ii).fvl=prctile(datout,50);
                    cntpt(ii).f05=prctile(datout,5);
                    cntpt(ii).f95=prctile(datout,95);
                    end
                end
                end
        end
        
        
        
%         vlsout{ii}=combvls(ptvlsind,:);
  
       
   end
end
end

%provide indices for bas and rec ind
for ii=1:length(avls.bastms)
   basbnds=datenum(avls.bastms{ii});
   wnbnds=datenum(avls.wntms{ii});
   recbnds=datenum(avls.rectms{ii});
   
   if(exist('dayout'))
    dvl.basind{ii}=find(dayout>=basbnds(1)&dayout<=basbnds(2));
    dvl.bas_sl(ii)=dvl.basind{ii}(end);
    dvl.wnind{ii}=find(floor(dayout)==(wnbnds(1)+avls.offset_dys(ii)));
    dvl.wn_sl(ii)=dvl.wnind{ii}(end);
    if(~isempty(recbnds)) 
        dvl.recind{ii}=find(dayout>=recbnds(1)&dayout<=recbnds(2));
        if(~isempty(dvl.recind{ii}))
            dvl.rec_sl(ii)=dvl.recind{ii}(end);
        else
            dvl.rec_sl(ii)=NaN;
        end
    else
        dvl.recind{ii}=[];
        dvl.rec_sl(ii)=NaN;
    end
   end
end



cmd=['cd ' avls.baspath ]
eval(cmd);
if(avls.DOCONTSYNANAL)
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat  cntseq vlsout outnotect dayout avls dvl'];

    eval(cmd);
    else
        cmd=['save sumdata.mat  cntseq vlsout outnotect dayout avls dvl'];
        eval(cmd)
    end
        
elseif(avls.DOCONTPTANAL)
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat cntpt'];

    eval(cmd);
    else
        cmd=['save sumdata.mat  cntpt'];
        eval(cmd)
    end
    
else
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat  outnotect dayout vlsout avls dvl'];
        eval(cmd); 
    else
        cmd=['save  sumdata.mat  outnotect dayout vlsout avls dvl'];
        eval(cmd); 
    end
end
    
end






