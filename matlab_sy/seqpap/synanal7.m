%CALL THIS FUNCTION WITH SPECIFIC INPUTS set in a bird-specific .m file
%Data saved to a .mat file
%synanal7 from synanal6: changed avls.ALLSEQNT to cell to allow different
%targets for different

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
    outallstruct.labels=labels(indnonzero);
    outallstruct.tms=alltms(indnonzero);
    
    %this is a hack to get inds for seqanalysis
    zerovec=zeros(length(labels),1);
    for ntind=1:length(avls.ALLSEQNT{crind})
       
        crnt=avls.NT{avls.ALLSEQNT{crind}(ntind)};
        ind=find(labels==crnt);
        zerovec(ind)=1;
    end 
    seqnz=find(zerovec);
    outseqstruct.labels=labels(seqnz);
    outseqstruct.tms=alltms(seqnz);
      cmd=['save  ' crbt '.mat outallstruct outseqstruct fvpt vls'];
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
combseqlabels=[];
combseqtms=[];
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
    comblabels=[comblabels;outallstruct.labels]
    combtms=[combtms;outallstruct.tms]
    
    combseqlabels=[combseqlabels;outseqstruct.labels]
    combseqtms=[combseqtms;outseqstruct.tms];
    ptvls=[];
    for crptnt=avls.ANALPTNT
        ptvls=[ptvls ;vls{crptnt}]
    end
    [out,ind]=sort(ptvls(:,1))
    combvls=[combvls; ptvls(ind,:)];
end


%now set up the individual responsevecs for syntax calculations.
%divide into daysed

    unq_days=unique(floor(combtms));
%looping through each exp. day;
if(~isempty(unq_days))
for ii=1:length(unq_days)
    if(ii==1)
        SEEDMEAN=0.5;
    end
   res_vec=[];
    crday=unq_days(ii);
   
       crdayind=find(floor(combtms)==crday)
       crseqdayind=find(floor(combseqtms)==crday)
       ptvlsind=find(floor(combvls(:,1))==crday)
   
    vlsout{ii}=combvls(ptvlsind,:);
   res_vec=zeros(length(crseqdayind),1);
   for ntind=1:length(avls.NT)
       
            indnt{ntind}=find(comblabels(crdayind)==avls.NT{ntind});
            outnotect(ii,ntind)=length(indnt{ntind})
            seqindnt{ntind}=find(combseqlabels(crseqdayind)==avls.NT{ntind})
       
           indnt{ntind}=find(comblabels(crdayind)==avls.NT{ntind});
            outnotect(ii,ntind)=length(indnt{ntind})
   end
  
   res_vec(seqindnt{avls.SEQTRGNT})=1;
   dayout(ii)=crday;
    
       
       if(avls.DOCONTSYNANAL||avls.DOCONTPTANAL)
        %rhis is the sequence analysis
        %res_vec is a series of zeros and ones; need to change expected
        %probability from 0.5 to prob. from prior day
        if(avls.DOCONTSYNANAL)
            if(~isempty(find(res_vec>0))&~isempty(find(res_vec~=1)))
                runanalysis(res_vec,1,SEEDMEAN);
                cmd=['load resultsindividual.mat'];
                eval(cmd);
            
                tmscomb=[combseqtms(crseqdayind)];
   %average values with the same "time", i.e. in same song    
            [cntseq(ii)]=ave_vls(tmscomb,p,p05,p95)
    
 
            SEEDMEAN=mean(cntseq(ii).vl);
            else
                SEEDMEAN=0.5;
            end
        end
        %this is the continuous pitch analysis.
        if(avls.DOCONTPTANAL)
            WIN=15;
            if(~isempty(vlsout{ii}))
                if(length(vlsout{ii}(:,1))>2*WIN)
%                     cntpt(ii).tms=runave2(vlsout{ii}(:,1),WIN);
                      startind=ceil(WIN/2);
                      tmptms=vlsout{ii}(startind:(end-ceil(WIN/2)),1)

                    [datout]=calcerbndscnt(vlsout{ii}(:,2),WIN);
                    if(~isempty(datout))
                        tmpvl=prctile(datout,50);
                        tmp05=prctile(datout,5);
                        tmp95=prctile(datout,95);
                        cntpt(ii)=ave_vls(tmptms,tmpvl,tmp05,tmp95);
                    end
                end
                
                
            end
           %average  
%            cntpt(ii).tms
%             
        end
        
        
        
%         vlsout{ii}=combvls(ptvlsind,:);
  
       
   end
end
end

%provide indices for bas and rec ind
for ii=1:length(avls.bastms)
   basbnds=datenum(avls.bastms{ii});
   wnbnds=datenum(avls.wntms{ii});
   if(~isempty(avls.rectms{ii}))
    recbnds=datenum(avls.rectms{ii});
   else
       recbnds=[];
   end
   if(exist('dayout'))
    dvl.basind{ii}=find(dayout>=basbnds(1)&dayout<=basbnds(2));
    dvl.bas_sl{ii}=dvl.basind{ii};
    dvl.wnind{ii}=find(floor(dayout)==(wnbnds(1)+avls.offset_dys(ii)));
    dvl.wn_sl(ii)=dvl.wnind{ii}(end);
    if(isfield(avls,'SPLITDAY'))
        if(avls.SPLITDAY)
            dvl.SPLITDAY=1;
        else
            dvl.SPLITDAY=0;
        
        end
    else
        dvl.SPLITDAY=0;
    end
    
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
        cmd=['save -append  sumdata.mat  cntseq vlsout outnotect combseqtms combseqlabels dayout avls dvl'];

    eval(cmd);
    else
        cmd=['save sumdata.mat  cntseq vlsout outnotect dayout avls dvl'];
        eval(cmd)
    end
end     
if(avls.DOCONTPTANAL)
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat cntpt'];

    eval(cmd);
    else
        cmd=['save sumdata.mat  cntpt'];
        eval(cmd)
    end
    
else
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat  outnotect dayout combseqtms combseqlabels vlsout avls dvl'];
        eval(cmd); 
    else
        cmd=['save  sumdata.mat  outnotect dayout vlsout combseqtms combseqlabels avls dvl'];
        eval(cmd); 
    end
end
else
    if(exist('sumdata.mat'))
        cmd=['save -append  sumdata.mat  avls '];
        eval(cmd); 
    else
        cmd=['save  sumdata.mat  avls '];
        eval(cmd); 
    end
end

function [datout]=ave_vls(tms,vls,vls05,vls95)

    unqtms=unique(tms);
    for tmind=1:length(unqtms)
                    selind=find(tms==unqtms(tmind));
                    datout.tms(tmind)=unqtms(tmind);
                    datout.vl(tmind)=mean(vls(selind));
                    datout.vl05(tmind)=mean(vls05(selind));
                    datout.vl95(tmind)=mean(vls95(selind));
   end


