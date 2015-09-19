%%assignvals_script

for ii=1:length(fbins)
    edges{ii}=[fbins{ii}(1):50:fbins{ii}(end)]
end


avls.usex=[0]

clear colvals
colvals{1}=muon;
colvals{2}=acon
colvals{3}=diron;

graphvals.numcol=3
graphvals.col='rkc'
% graphvals.plttext=1


avls.tmon=tmn;
avls.tmoff=tmf;
graphvals.tmn=tmn;
graphvals.tmf=tmf;
graphvals.acon=acon;

avls.edges=edges;

graphvals.colvals=colvals
graphvals.chunkdata=1

if exist('prefix')
    [pathvl]=addprefix(pathvl,prefix);
end
avls.pvls=pathvl
avls.datfile=matfilename;
avls.cvl=catchvl
avls.sumpath=sumpath
avls.mtflnm=matfilename  
avls.NT=NT
avls.acon=acon;
avls.muon=muon;
avls.NFFT=NFFT
avls.fbins=fbins
avls.tshft=tbinshft
avls.numnotes=numnotes
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.deadtm=.16/24;
avls.muoffset=2/24;
avls.acoffset=1/24;
avls.pretm=5/24

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


