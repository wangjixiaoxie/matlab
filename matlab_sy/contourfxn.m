%8.23.09
%originally created to take a ptst (pitchstruct),
%called by bird-specific scripts(i.e. contouranalbk59w37.m

function [] = contourfxn(ptst)
avls=ptst.avls;

for ii=1:length(ptst.muinds)
    crind=ptst.muinds(ii); 
    pathvl=avls.pvls{crind}
    cmd=['cd ' pathvl];eval(cmd);
    bt=avls.cvl{crind};
    NT=avls.NT{crind};PRENT='';PSTNT='';
    fvpt=findwnote4(bt,NT,PRENT,PSTNT,ptst.pt_tbinshft,ptst.fbins,ptst.pt_NFFT,1,'obs0');
    fv=findwnote4(bt,NT,PRENT,PSTNT,ptst.con_tbinshft,ptst.conbins,ptst.con_NFFT,1,'obs0');
    pitchdata{crind}=jc_pitchcontourFV(fv,1024,1020,1,ptst.conbins(1), ptst.conbins(2), 3, 'obs0');
    contours=pitchdata{crind};
    tms=maketimevec2(bt);
    initsamptime=512/32000;
    initsamptimediff=1/8000;
    pitchtms=initsamptime:initsamptimediff:(4096-512)/32000
    pitchtms=pitchtms+ptst.con_tbinshft
    cmd=['save ' bt '.mat contours tms fvpt pitchtms'];
    eval(cmd);
end
