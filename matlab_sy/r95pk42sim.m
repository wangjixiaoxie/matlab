%use X.tmp file with final template, to simulate hits and misses.


pathvl{1}='/doyale/twarren/r95pk42/templteste2/';
catchvl{1}='batch1920comb'
pathvl{2}='/doyale/twarren/r95pk42/templteste3/'
catchvl{2}='batch.catch.keep'
pathvl{3}='/doya/r95pk42/wnon2/'
catchvl{3}='batch2223'
pathvl{4}='/doya/r95pk42/wnon2/'
catchvl{4}='batch24.catch.keep'
pathvl{5}=pathvl{4}
catchvl{5}='batch25.catch.keep'
pathvl{6}='/doyale/twarren/r95pk42/wnon3/'
catchvl{6}='batch26.keep.catch'
pathvl{7}='/doyale/twarren/r95pk42/wnon3/'
catchvl{7}='batch27.keep.catch.rand'

valscomb=[];
fvcomb=[]
for ii=1:length(pathvl)
    strcmd=['cd '  pathvl{ii}];
    eval(strcmd);
    bt=catchvl{ii}
    
    tbinshft=0.005;
    NFFT=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    fbins=[2000,5000];
    save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
   load BINS_B
    NT='e';PRENT='d';PSTNT='';

    
       % edges=[6000:75:8000];
    
    fv=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0',1);
    fvcomb=[fvcomb fv]
    vals=getvals(fv,1,'TRIG');

    valscomb=[valscomb;vals]
end