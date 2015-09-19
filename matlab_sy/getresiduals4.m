%this function is used for short pulse analysis across birds.
%START OFF BY SETTING pitchtms, assuming they are constant across birds.
function [os,combfvst,combstcr]= getresiduals(sps,sts,ind) 

RESID_TMS=[.070 .1];
OFF_TMS=[0.078 0.088];
TMBUFF=0.1
os=struct([]);

for spnum=1:length(ind)
    spsvl=ind(spnum)
    crsp=sps(spsvl);
    for grp=1:length(crsp.datbnds)
        %generate mean catch contour
        
        mnct=getmeancatch(sts, crsp, grp)
        
        stbnds=crsp.datbnds{grp}
        %generate a matrix of residuals
        [os,combfvst,combstcr]=getresid(sts,spnum, crsp,grp,mnct, RESID_TMS,OFF_TMS,stbnds,TMBUFF,os)
    end
end
        
function [mnct] = getmeancatch(sts,crsp,grp)
    stind=crsp.stind;
    path=sts(stind).path
    mfile=sts(stind).matfilename;
    pvlctind=crsp.catch{grp}
    [avls]=load_avls(path,mfile)
    crctcomb=[];
    for ctind=1:length(pvlctind)
        crct=pvlctind(ctind);
        pvlct=crsp.inds(crct);
        [contours,crctind, crfbind,fvst,pitchtms]=loadcontours(avls, pvlct);
        crctcontours=contours(:,crctind);
        crctcomb=[crctcomb crctcontours];
    end
    mnct=mean(crctcomb,2);
    
function [os,combfvst,combstcr] = getresid(sts,spnum,crsp,grp,mnct,RESID_TMS,OFF_TMS,stbnds,TMBUFF,os)
    stind=crsp.stind;
    path=sts(stind).path
    mfile=sts(stind).matfilename;
    pvlfbind=crsp.fb{grp}
    pvlfb=crsp.inds(pvlfbind);
    [avls]=load_avls(path,mfile);
    [combfb,combfvst,pitchtms,del_list,TMDIFF,combstcr,combctr,combct]=combinefbcontours(avls, pvlfb);
    [residmatfb,residmatfboff,residmatct,residmatctoff,residtmsct,residtmsctoff,residtmsfb,residtmsfboff,combctcon,combfbcon]=mkresiduals(pitchtms,combfb,combfvst,TMBUFF,RESID_TMS,OFF_TMS,mnct,del_list,stbnds,TMDIFF,combctr,combct);
    [tm,mndiff,mnctlim,mner,mnresid,steresid]=calcoff(mnct,residmatfboff,residmatctoff,residtmsfboff,TMDIFF);
    os(spnum).outmatfb{grp}=residmatfb;
    os(spnum).outmatct{grp}=residmatct
    os(spnum).ct_tms{grp}=residtmsct;
    os(spnum).fb_tms{grp}=residtmsfb;
    os(spnum).tm{grp}=tm;
    os(spnum).mndiff{grp}=mndiff;
    os(spnum).mnct{grp}=mnctlim;
    os(spnum).mner{grp}=mner;
    os(spnum).mnresid{grp}=mnresid;
    os(spnum).steresid{grp}=steresid;
    os(spnum).steresid{grp}=steresid;
    os(spnum).combctcon{grp}=combctcon
    os(spnum).combfbcon{grp}=combfbcon;
    
    
    
    function [avls]=load_avls(path,mfile)
    cmd=['cd ' path];
    eval(cmd);
    cmd=['load ' mfile];
    eval(cmd);

function [contours,crctind,crfbind,fvst,pitchtms]=loadcontours(avls,crind)
    if(isfield(avls,'baspath'))
        path=[avls.baspath avls.pvls{crind}]
    else
        path=[avls.pvls{crind}];
    end
    cmd=['cd ' path];
    eval(cmd);
    bt=avls.cvl{crind}
    cmd=['load ' bt '.mat']
    eval(cmd);
    
function [combfb,combfvst,pitchtms, del_list,TMDIFF,combstcr,combctr,combct]=combinefbcontours(avls,pvlfbind)
    combfvst=[]
    combfb=[];
    combstcr=[];
    combctr=[];
    del_list=[];
    for fbind=1:length(pvlfbind)
        crfb=pvlfbind(fbind);
        [contours,crctind,crfbind,fvst,pitchtms]=loadcontours(avls,crfb);
        crdel=avls.del{crfb}*ones(length(fvst),1);
        
        del_list=[del_list crdel'];
        
        combfvst=[combfvst fvst];
        fvcr=fvst;
        for ii=1:length(fvcr)
            stcr(ii)=fvcr(ii).STIMTIME
        end
        combstcr=[combstcr stcr/1000+avls.del{crfb}]
        combfb=[combfb contours(:,crfbind)];
        combct=[combfb contours(:,crctind)];
        combctr=[combctr contours];
    end
    pitchtms=.066:1.25e-4:.162
    TMDIFF=pitchtms(2)-pitchtms(1)
    

%input are a matrix of feedback contours, a matric of matching fvst
%structs, pitchtms, TMBUFF, and the mean catch contour
%this code is borrowed from stimvardelay.m
    
function [residmatfb,residmatfboff,residmatct,residmatctoff,mat_tmsct,mat_tmsctoff,mat_tmsfb,mat_tmsfboff,combctcon,combfbcon]=mkresiduals(pitchtms,combfb,combfvst,TMBUFF,RESID_TMS,OFF_TMS,mnct, del_list,stbnds,TMDIFF,combctr,combct)            
    combctcon=[];
    combfbcon=[];
    fbind=[];
    ctind=[];
    tmbnds=get_tmbnds(pitchtms,RESID_TMS);
    tmbndsoff=get_tmbnds(pitchtms,OFF_TMS);
    mnct=mnct(tmbnds(1):tmbnds(2));
    mnctoff=mnct(tmbndsoff(1):tmbndsoff(2));
    %each ROW is a different diff contour??
    for ii=1:length(combctr(1,:))
        diff_fbmat(:,ii)=combctr(tmbnds(1):tmbnds(2), ii)-mnct;
        diff_fbmatoff(:,ii)=combctr(tmbndsoff(1):tmbndsoff(2),ii)-mnctoff;
    end
    [stcr]=getstimtimes(combfvst);
    actual_st=stcr/1000+del_list;
    stinds=find(actual_st>=stbnds(1)&actual_st<=stbnds(2));
    
    for jj=1:length(stinds)
        crind=stinds(jj);
        if(combfvst(crind).STIMCATCH==1)
            [combctcon]=[combctcon combctr(:,crind)] 
            ctind=[ctind crind];
        elseif(combfvst(crind).STIMCATCH==0)
            [combfbcon]=[combfbcon combctr(:,crind)];
            fbind=[fbind crind];
        end
   end    
    %these should be fb stim times
    crfbmat=diff_fbmat(:,fbind);
    crfboffmat=diff_fbmatoff(:,fbind);
    
    crctmat=diff_fbmat(:,ctind);
    crctoffmat=diff_fbmatoff(:,ctind);
    residfb_tms=actual_st(fbind)
    residct_tms=actual_st(ctind);
    
    [residmatfb,residmatfboff,residmatct,residmatctoff,mat_tmsct,mat_tmsctoff,mat_tmsfb,mat_tmsfboff]=mkadjmat(crfbmat,crfboffmat,crctmat,crctoffmat,residfb_tms,residct_tms,tmbnds,RESID_TMS,OFF_TMS,TMDIFF,TMBUFF);
    
function [residmatfb,residmatfboff,residmatct,residmatctoff,mat_tmsct,mat_tmsctoff,mat_tmsfb,mat_tmsfboff]=mkadjmat(crfbmat,crfboffmat,crctmat,crctoffmat,residfb_tms,residct_tms,tmbnds,RESID_TMS,OFF_TMS,TMDIFF,TMBUFF,fbind,ctind)
    residmatfb=zeros(1000,length(crfbmat(1,:)));
    residmatfboff=zeros(1000,length(crfboffmat(1,:)));
    residmatct=zeros(1000,length(crctmat(1,:)));
    residmatctoff=zeros(1000,length(crctoffmat(1,:)));
    lnmat=length(crfbmat(:,1));
    lnmatoff=length(crfboffmat(:,1));
    
    for ii=1:length(residfb_tms)
        adjtmshft(ii)=RESID_TMS(1)-residfb_tms(ii);
        ptshft(ii)=round(adjtmshft(ii)/TMDIFF);
        ptbuff=TMBUFF/TMDIFF
        residmatfb(floor(ptshft(ii)+ptbuff):floor(ptbuff+ptshft(ii)+lnmat-1),ii)=crfbmat(:,ii);
      
        adjoff(ii)=OFF_TMS(1)-residfb_tms(ii);
        ptshftoff(ii)=round(adjoff(ii)/TMDIFF);
        residmatfboff(floor(ptshftoff(ii)+ptbuff):floor(ptbuff+ptshftoff(ii)+lnmatoff-1),ii)=crfboffmat(:,ii);
        
        
    end

    for ii=1:length(residct_tms)
        adjtmshft(ii)=RESID_TMS(1)-residct_tms(ii);
        ptshft(ii)=round(adjtmshft(ii)/TMDIFF);
        ptbuff=TMBUFF/TMDIFF
        residmatct(floor(ptshft(ii)+ptbuff):floor(ptbuff+ptshft(ii)+lnmat-1),ii)=crctmat(:,ii);
         
        adjoff(ii)=OFF_TMS(1)-residct_tms(ii);
        ptshftoff(ii)=round(adjoff(ii)/TMDIFF);
        residmatctoff(floor(ptshftoff(ii)+ptbuff):floor(ptbuff+ptshftoff(ii)+lnmatoff-1),ii)=crctoffmat(:,ii);
        
        
    end
    mat_tmsct=-TMBUFF:TMDIFF:TMDIFF*length(residmatct(:,1))-TMBUFF;
    mat_tmsct=mat_tmsct(1:end-1)
    
    mat_tmsctoff=-TMBUFF:TMDIFF:TMDIFF*length(residmatctoff(:,1))-TMBUFF;
    mat_tmsctoff=mat_tmsctoff(1:end-1)
    
    mat_tmsfb=-TMBUFF:TMDIFF:TMDIFF*length(residmatfb(:,1))-TMBUFF;
    mat_tmsfb=mat_tmsfb(1:end-1)
    
    mat_tmsfboff=-TMBUFF:TMDIFF:TMDIFF*length(residmatfboff(:,1))-TMBUFF;
    mat_tmsfboff=mat_tmsfboff(1:end-1)

function [tmbnds]=get_tmbnds(pitchtms,PT_TMS)
    for ii=1:2
        difftms=pitchtms-PT_TMS(ii);
        minvl=min(abs(difftms));
        indmin=find(abs(difftms)==minvl)
        tmbnds(ii)=indmin;
    end
    
    function [stcr]=getstimtimes(combfvst,del_list)
        for ii=1:length(combfvst)
            stcr(ii)=combfvst(ii).STIMTIME
        end

%this function rewritten 11.9.09 to only measure in time window set by
%OFST_TMS

function   [tm,meanout,mnctlim,meaner,meanresid,steresid,residtmsfb ]=calcoff(mnct,residmatfb,residmatct,residtmsfb,TMDIFF)
        
        %go through residmat and find residmat~=0
        lnfb=length(residmatfb(:,1));
        lnct=length(residmatct(:,1));
        numfb=length(residmatfb(1,:));
        if(lnct<lnfb)
            ln=lnct
%             tms=ct_tms
        else
            ln=lnfb
%             tms=fbtms
        end
            
        for ii=1:ln
            nzindfb{ii}=find(residmatfb(ii,:)~=0);
            nzindct{ii}=find(residmatct(ii,:)~=0);
            if(~isempty(nzindfb))
                nonzeroind(ii)=length(nzindfb{ii});
            else
                nonzeroind(ii)=0;
            end
        end
        suffind=find(nonzeroind>.5*numfb);
        startind=suffind(1);
        tm=residtmsfb(startind);
        endind=startind+.008/TMDIFF;
        inds=startind:endind;
        
        for ii=1:length(inds)
            crind=inds(ii);
            meanresid(ii)=mean(residmatfb(crind,nzindfb{crind}),2)-mean(residmatct(crind,nzindct{crind}),2);
            steresid(ii)=std(residmatfb(crind,nzindfb{crind}),0,2)./sqrt(length(nzindct{crind}));
        end
            meanout=mean(meanresid);
            meaner=mean(steresid);
            mnctlim=mean(mnct(1:65));
        %calculate meanbas
        

        
        