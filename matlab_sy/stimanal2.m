function [avls] = stimanal(avls)
muinds=[]
clear valsa
SOUNDSTR='wav'
STIMSTR='rig'

REPSOUNDSTR='und'
REPSTIMSTR='im'

REPINDS=[]

%check if data directory exists; if it doesn't, create it.
if(isfield(avls,'baspath'))
    cmd=['cd ' avls.baspath]
    eval(cmd);
end
if (exist('avls.datdir')~=7)
    cmd=['mkdir ' avls.datdir]
    eval(cmd);
end


for ii=1:length(avls.ptfv)
    crind=avls.ptfv(ii);
    pathvl=avls.pvls{crind}
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]    
    else
        cmd=['cd ' pathvl]
    end
    eval(cmd);
    bt=avls.cvl{crind};
    cmd=['load ' bt '.mat']
    eval(cmd);
        tbinshft=-.02;
 
    NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
     if(iscell(avls.fbins))
        crfbins=avls.fbins{crind}
    else
        crfbins=avls.fbins;
    end

    % frequency analysis just for 'b'
    if(isfield(avls,'CORSOUNDSTR'))
        if(~isempty(intersect(avls.CORSTRINDS,crind)))
            SOUNDSTR=avls.CORSOUNDSTR
            STIMSTR=avls.CORSTIMSTR
        else
            SOUNDSTR=avls.SOUNDSTR
            STIMSTR=avls.STIMSTR;
        end
    else
            SOUNDSTR=avls.SOUNDSTR
            STIMSTR=avls.STIMSTR;
    end
    
    
    NT=avls.NT{crind};PRENT='';PSTNT='';
    fvst=findwnote9(bt,NT,PRENT,PSTNT,tbinshft,crfbins,NFFT,1,'obs0',0,avls.STIMBNDS, SOUNDSTR,STIMSTR);
% %     fv=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0',0)
    vals=getvals2(fvpt,1,'TRIG');
    
    cmd=['save -append ' bt '.mat vals fvst'];
    eval(cmd);
end



clear notind fbind ctind
for ii=1:length(avls.catchstimfv)
run_num=avls.catchstimfv(ii);
   notind{run_num}=[]
    fbind{run_num}=[];
    ctind{run_num}=[];
    
    pathvl=avls.pvls{run_num}
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end
    eval(cmd);
    bt=avls.cvl{run_num}
    cmd=['load ' bt '.mat']  ;
    eval(cmd);
    for ii=1:length(fvst)
        if(fvst(ii).STIMTRIG)
            if(fvst(ii).STIMCATCH)
                ctind{run_num}=[ctind{run_num} ii]
            else
                fbind{run_num}=[fbind{run_num} ii]
            end
        else
            notind{run_num}=[notind{run_num} ii]
        end
    end
    crfbind=fbind{run_num};
    crnotind=notind{run_num};
    crctind=ctind{run_num};
    pathvl=avls.pvls{run_num}
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end
   

    eval(cmd);
   
    %these values are useful for plotting.
    edges=avls.HST_EDGES
    ctvls=vals(crctind,2);
    fbvls=vals(crfbind,2);
    ctmean=mean(ctvls);
    fbmean=mean(fbvls);
    stdfb=std(fbvls);
    stdct=std(ctvls);
    hstoutctind=histc(ctvls,edges);
    hsctnrm=hstoutctind./length(crctind);
    hstoutfbind=histc(fbvls,edges);
    hsfbnrm=hstoutfbind./length(crfbind);
   
    clear stct
   
    cmd=['save -append ' bt '.mat crfbind crnotind crctind hsctnrm hsfbnrm  ctmean fbmean stdfb stdct vals ']  ;
    eval(cmd);
end

for ii=1:length(avls.analfv)
    run_num=avls.analfv(ii);
    pathvl=avls.pvls{run_num}
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd '  pathvl]
    end
        eval(cmd)
    cmd=['load ' avls.cvl{run_num} '.mat'];
    eval(cmd);
    
    clear stcr
    for jj=1:length(fvst)
       stcr(jj)=fvst(jj).STIMTIME;
    end
    [avls.adjtimes(run_num,1), avls.adjtimes(run_num,2)]=get_times(fvst);
    
    if(isfield(avls,'mnbas'))
        [avls.catchz(run_num),avls.stimz(run_num)]=mk_zvls(ctmean, fbmean, avls);
    end
    if(isfield(avls,'wnin'))
        [avls.wn,avls.wnrev]=getwntimes(avls.wnin, avls.wnrevin);
    end
    
    avls.stcr{run_num}=stcr;
    avls.ctmean(run_num)=ctmean
    avls.fbmean(run_num)=fbmean
    avls.stdfb(run_num)=stdfb
    avls.stdct(run_num)=stdct
    ptvls=getvals_sec(fvpt,1,'trig')
    avls.ptvls{run_num}=ptvls;
    avls.hsctnrm{run_num}=hsctnrm
    avls.hsfbnrm{run_num}=hsfbnrm
    avls.crctind{run_num}=crctind
    avls.crfbind{run_num}=crfbind
    if(isfield(avls,'REMOVEOUTLIERS'))
       if(avls.REMOVEOUTLIERS)
           crfbindlim=remove_outliers(avls,run_num);
           avls.crfbindlim{run_num}=crfbindlim;
           fbvls=vals(crfbindlim,2);
           avls.fbmeanlim(run_num)=mean(fbvls);
           avls.stdfblim(run_num)=std(vals(crfbindlim,2));
           
           avls.hsfbnrm{run_num}=histc(fbvls,avls.HST_EDGES)./length(fbvls);
       end
    end
    cmd=['save -append ' bt '.mat crfbindlim']  ;
end
    
if(isfield(avls,'baspath'))
    cmd=['cd ' avls.baspath avls.datdir]
else
    cmd=['cd ' avls.datdir]
end




eval(cmd);
cmd=['save datsum.mat avls']
eval(cmd)

function [crfblim]=remove_outliers(avls,run_num)
    NUMSTD=avls.NUMSTD;
    
    crfbind=avls.crfbind{run_num};
    stcr=avls.stcr{run_num};
   
    mnvl=mean(stcr);
    stdvl=std(stcr);
    minst=mnvl-NUMSTD*stdvl;
    maxst=mnvl+NUMSTD*stdvl;
    stlim=find(stcr>minst&stcr<maxst);
    crfblim=intersect(crfbind,stlim);
   


function [ctz,stimz] = mk_zvls(ctmean, fbmean,avls)
    ctz=(ctmean-avls.mnbas)/avls.stdbas;
    stimz=(fbmean-avls.mnbas)/avls.stdbas;
    
function [tminit, tmfinal]=get_times(fvst)
    inda=find(fvst(1).fn=='_')
    indb=find(fvst(1).fn=='.')
    if(indb(1)-inda(2)~=7)
        tminit=fn2datenum(fvst(1).fn);
        tmfinal=fn2datenum(fvst(end).fn);
    else
        tminit=fn2datenumsec(fvst(1).fn);
        tmfinal=fn2datenumsec(fvst(end).fn);
    end