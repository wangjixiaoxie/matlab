%this function is used for short pulse analysis across birds.
%START OFF BY SETTING pitchtms, assuming they are constant across birds.
function [outstr,combfvst,combstcr]= getresiduals(sps,sts,ind) 

PT_TMS=[.078 .100]
TMBUFF=0.1

for spsnum=1:length(ind)
    spsvl=ind(spsnum)
    crsp=sps(spsvl);
    for grp=1:length(crsp.datbnds)
        %generate mean catch contour
        
        mnct=getmeancatch(sts, crsp, grp)
        
        stbnds=crsp.datbnds{grp}
        %generate a matrix of residuals
        [outstr(spsnum).outmat{grp}, outstr(spsnum).out_tms{grp},outstr(spsnum).tm{grp},outstr(spsnum).mndiff{grp},outstr(spsnum).mnct{grp},outstr(spsnum).mner{grp},outstr(spsnum).mnresid{grp},outstr(spsnum).steresid{grp},combfvst,combstcr]=getresid(sts,crsp, grp,mnct, PT_TMS,stbnds,TMBUFF)
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
    
function [residmat,residtms,tm,mndiff,mnctlim,mner,mnresid,steresid,combfvst,combstcr] = getresid(sts,crsp,grp,mnct,PT_TMS,stbnds,TMBUFF)
    stind=crsp.stind;
    path=sts(stind).path
    mfile=sts(stind).matfilename;
    pvlfbind=crsp.fb{grp}
    pvlfb=crsp.inds(pvlfbind);
    [avls]=load_avls(path,mfile);
    [combfb,combfvst,pitchtms,del_list,TMDIFF,combstcr]=combinefbcontours(avls, pvlfb);
    [residmat,residtms]=mkresiduals(pitchtms,combfb,combfvst,TMBUFF,PT_TMS,mnct,del_list,stbnds,TMDIFF);
    [tm,mndiff,mnctlim,mner,mnresid,steresid]=calcoff(mnct,residmat,residtms,crsp,TMDIFF);
    
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
    
function [combfb,combfvst,pitchtms, del_list,TMDIFF,combstcr]=combinefbcontours(avls,pvlfbind)
    combfvst=[]
    combfb=[];
    combstcr=[];
    del_list=[];
    for fbind=1:length(pvlfbind)
        crfb=pvlfbind(fbind);
        [contours,crctind,crfbind,fvst,pitchtms]=loadcontours(avls,crfb);
        crdel=avls.del{crfb}*ones(length(fvst(crfbind)),1);
        del_list=[del_list crdel'];
        
        combfvst=[combfvst fvst(crfbind)];
        fvcr=fvst;
        for ii=1:length(fvcr)
            stcr(ii)=fvcr(ii).STIMTIME
        end
        combstcr=[combstcr stcr/1000+avls.del{crfb}]
        combfb=[combfb contours(:,crfbind)];
    end
    pitchtms=.066:1.25e-4:.162
    TMDIFF=pitchtms(2)-pitchtms(1)
    

%input are a matrix of feedback contours, a matric of matching fvst
%structs, pitchtms, TMBUFF, and the mean catch contour
%this code is borrowed from stimvardelay.m
    
function [residmat,mat_tms]=mkresiduals(pitchtms,combfb,combfvst,TMBUFF,PT_TMS,mnct, del_list,stbnds,TMDIFF)            
    tmbnds=get_tmbnds(pitchtms,PT_TMS);
    mnct=mnct(tmbnds(1):tmbnds(2));
    %each ROW is a different diff contour??
    for ii=1:length(combfb(1,:))
        diff_fbmat(:,ii)=combfb(tmbnds(1):tmbnds(2), ii)-mnct;
    end
    [stcr]=getstimtimes(combfvst);
    actual_st=stcr/1000+del_list;
    stinds=find(actual_st>=stbnds(1)&actual_st<=stbnds(2));
    crfbmat=diff_fbmat(:,stinds);
    crstimes=actual_st(stinds);
    
    [residmat,mat_tms]=mkadjmat(crfbmat,crstimes,tmbnds,PT_TMS,TMDIFF,TMBUFF);
    
function [residmat,mat_tms]=mkadjmat(crfbmat,residtms,tmbnds,PT_TMS,TMDIFF,TMBUFF)
    residmat=zeros(1000,length(crfbmat(1,:)));
    lnmat=length(crfbmat(:,1));
   
    for ii=1:length(residtms)
        adjtmshft(ii)=PT_TMS(1)-residtms(ii);
        ptshft(ii)=round(adjtmshft(ii)/TMDIFF);
        ptbuff=TMBUFF/TMDIFF
        residmat(floor(ptshft(ii)+ptbuff):floor(ptbuff+ptshft(ii)+lnmat-1),ii)=crfbmat(:,ii);
    end
    mat_tms=-TMBUFF:TMDIFF:TMDIFF*length(residmat(:,1))-TMBUFF;
    mat_tms=mat_tms(1:end-1)

function [tmbnds]=get_tmbnds(pitchtms,PT_TMS)
    for ii=1:2
        difftms=pitchtms-PT_TMS(ii);
        minvl=min(abs(difftms));
        indmin=find(difftms==minvl)
        tmbnds(ii)=indmin;
    end
    
    function [stcr]=getstimtimes(combfvst,del_list)
        for ii=1:length(combfvst)
            stcr(ii)=combfvst(ii).STIMTIME
        end
    function   [tm,meanout,mnctlim,meaner,meanresid,steresid,residtms ]=calcoff(mnct,residmat,residtms,crsp,TMDIFF)
        
        %go through residmat and find residmat~=0
        ln=length(residmat(1,:));
        for ii=1:length(residmat(:,1))
            nzind{ii}=find(residmat(ii,:)~=0);
            if(~isempty(nzind))
                nonzeroind(ii)=length(nzind{ii});
            else
                nonzeroind(ii)=0;
            end
        end
        suffind=find(nonzeroind>.5*ln);
        startind=suffind(1);
        tm=residtms(startind);
        endind=startind+.008/TMDIFF;
        inds=startind:endind;
        
        for ii=1:length(inds)
            crind=inds(ii);
            meanresid(ii)=mean(residmat(crind,nzind{crind}),2);
            steresid(ii)=std(residmat(crind,nzind{crind}),0,2)./sqrt(length(nzind{crind}));
        end
            meanout=mean(meanresid);
            meaner=mean(steresid);
            mnctlim=mean(mnct(1:65));
        %calculate meanbas
        

        
        