% May 2009
% classification of D vs. UD
% Mimi Kao's ZF data
    % Only three birds: 1-2, 3-6, 7-10
    
% Step 1: find middle of note
    for j=1:10
        long=size(MKZF(j).pitchPRE,1);
        steps=round((long-100)/100);
        cvc=std(MKZF(j).pitchPRE');
        localCV=[];
        for i=1:82
            localCV(i)=mean(cvc(steps*i:steps*i+250));
        end
        [a,ind]=min(localCV);
        start(j)=ind*steps;
        figure;plot(MKZF(j).pitchPRE(start(j):start(j)+250,:))
    end
    start(6)=480;
% Step 2: make sure Dir song reduces CV - they all do
    for i=1:10
        CVratio(i)=std(MKZF(i).pitchDIR(start(i):start(i)+250,:)')/std(MKZF(i).pitchPRE(start(i):start(i)+250,:)');
    end
% Step 3: get residuals
    for i=1:10
        Resid(i).Dir=jc_residuals(MKZF(i).pitchDIR(start(i):start(i)+250,:));
        Resid(i).UnDir=jc_residuals(MKZF(i).pitchPRE(start(i):start(i)+250,:));
    end
% Step 4: put data into matrix for PCA
    Alldata=[];
    for i=1:10
        Alldata=[Alldata Resid(i).UnDir];
    end
    UDend=size(Alldata,2);
    for i=1:10
        Alldata=[Alldata Resid(i).Dir];
    end
% Step 5: do PCA
    [PC,SCORE,LATENT]=princomp(Alldata');
    [bestmin,marker,combo3]=PCAclassify2(SCORE,UDend);
% PCs look like sinusoids of different periods
% vco them!
% 1. take the raw data over the middle
i=1;
rawchunk=MKZF(i).rawDIR(:,(start(i)*4+512):(start(i)+250)*4+512);
 rawchunkPRE=MKZF(i).rawPRE(:,(start(i)*4+512):(start(i)+250)*4+512);
    