%%% 8.20.09
%%% MSB controls and additional analyses
% Zebra Finch notes
    % w5r22 - cardinal5 - late March 2009
    % exclusion of outliers
        ZFdata(1).Learned=mean(pitchshifted(:,[1:25 27:68 70:71 73:80 82:end])')-mean(pitch428(:,[1:40 42:end])');
        %%
                    for j=1
                        j
                        if isequal('up',ZFdata(j).direction);value=60;else value=40;end
                            CSsZF(j).data=ContingSim2(ZFdata(j).Targeting,ZFdata(1).ResidAC,value);
                        end
                    end
    % g5o78
    % December
    % stationary mean w/o outliers
    % cardinal/g5o78/Alldata1210.mat
        ZFdata(2).ResidAC=jc_residuals(pitchBase(1:1870,[80:181 183 185:220]));

        %
        pitchShift=[pitch1207A(1:1870,:),pitch1207B(1:1870,:),pitch1207C(1:1870,:),pitch1208A(1:1870,:),pitch1208B(1:1870,:),pitch1208C(1:1870,:),pitch1209A(1:1870,:),pitch1209B(1:1870,:)];
        ind=find((std(pitchShift(550:800,:)))<110)
        pitchShift=pitchShift(:,ind);
        ind2=find((std(pitchBase(440:1000,:)))<70)
        pitchBase=pitchBase(:,ind2);
        
%%% Discredit Model 1
        counter=zeros(1,6);
        for i=1:28
            edgeleft=round(median(Predict(i).Targeting)-2*std(Predict(i).Targeting));
            edgeright=round(median(Predict(i).Targeting)+2*std(Predict(i).Targeting));
            if edgeleft>Predict(i).onset
                counter(1)=counter(1)+1;
                neg0(counter(1))=Predict(i).LearnedNorm(edgeleft);
                if edgeleft-80>Predict(i).onset
                    counter(2)=counter(2)+1;
                    neg10(counter(2))=Predict(i).LearnedNorm(edgeleft-80);
                    if edgeleft-160>Predict(i).onset
                        counter(3)=counter(3)+1;
                        neg20(counter(3))=Predict(i).LearnedNorm(edgeleft-160);
                    end
                end

            end
            if edgeright<Predict(i).offset
                counter(4)=counter(4)+1;
                pos0(counter(4))=Predict(i).LearnedNorm(edgeright);
                if edgeright+80<Predict(i).offset
                    counter(5)=counter(5)+1;
                    pos10(counter(5))=Predict(i).LearnedNorm(edgeright+80);
                    if edgeright+160<Predict(i).offset
                        counter(6)=counter(6)+1;
                        pos20(counter(6))=Predict(i).LearnedNorm(edgeright+160);
                    end
                end

            end

        end
        counter=zeros(1,6);
        for i=30:31
            edgeleft=round(median(Predict(i).Targeting)-2*std(Predict(i).Targeting));
            edgeright=round(median(Predict(i).Targeting)+2*std(Predict(i).Targeting));
            if edgeleft>Predict(i).onset
                counter(1)=counter(1)+1;
                ZFneg0(counter(1))=Predict(i).LearnedNorm(edgeleft);
                if edgeleft-80>Predict(i).onset
                    counter(2)=counter(2)+1;
                    ZFneg10(counter(2))=Predict(i).LearnedNorm(edgeleft-80);
                    if edgeleft-160>Predict(i).onset
                        counter(3)=counter(3)+1;
                        ZFneg20(counter(3))=Predict(i).LearnedNorm(edgeleft-160);
                    end
                end

            end
            if edgeright<Predict(i).offset
                counter(4)=counter(4)+1;
                ZFpos0(counter(4))=Predict(i).LearnedNorm(edgeright);
                if edgeright+80<Predict(i).offset
                    counter(5)=counter(5)+1;
                    ZFpos10(counter(5))=Predict(i).LearnedNorm(edgeright+80);
                    if edgeright+160<Predict(i).offset
                        counter(6)=counter(6)+1;
                        ZFpos20(counter(6))=Predict(i).LearnedNorm(edgeright+160);
                    end
                end

            end

        end
        
        figure;plot(0,1,'.','MarkerSize',20,'Color','b')
        hold on;plot(-3,neg20,'.','MarkerSize',20,'Color','b')
        plot(-2,neg10,'.','MarkerSize',20,'Color','b')
        plot(-1,neg0,'.','MarkerSize',20,'Color','b')
        plot(1,pos0,'.','MarkerSize',20,'Color','b')
        plot(2,pos10,'.','MarkerSize',20,'Color','b')
        plot(3,pos20,'.','MarkerSize',20,'Color','b')
        % incorporate ZF experiments
        plot(-3,ZFneg20,'.','MarkerSize',20,'Color','r')
        plot(-2,ZFneg10,'.','MarkerSize',20,'Color','r')
        plot(-1,ZFneg0,'.','MarkerSize',20,'Color','r')
        plot(1,ZFpos0,'.','MarkerSize',20,'Color','r')
        plot(2,ZFpos10,'.','MarkerSize',20,'Color','r')
        plot(3,ZFpos20,'.','MarkerSize',20,'Color','r')

        plot([-3 -2 -1 0 1 2 3],[mean([neg20 ZFneg20]) mean([neg10 ZFneg10]) mean([neg0 ZFneg0]) 1 mean([pos0 ZFpos0]) mean([pos10 ZFpos10]) mean([pos20 ZFpos20])],'k')

% Next syllable
    % w5r22
        % base=428
        % shift=tmp2ampon509
        % B is not labeled but is 75ms after end of A
    % pk37bk19
        %base=1015files
        %shift=1021files
        % B  is labeled and follows 100ms after end of A and is not very stacky, but kind
        % of stacky
    % bk50w18
        % base =1005Afiles
        % shift=
%%%%% SUMMARY PLOT FOR SURROUNDING SYLLABLE DATA
    % Data in cardinal/Adjsylls.mat
    figure;hold on;

        for i=1:12
            aac(i)=Adjsylls(i).ccorr;
            bac(i)=Adjsylls(i).Bshift/Adjsylls(i).Ashift;
        end
        corrcoef(aac,bac)
        figure;hold on;plot(aac,bac,'.','MarkerSize',20)
        hold on;plot([-0.3 0.7],[-0.3 0.7])
        hold on;plot([0 0],[-0.3 0.7],'k')
        hold on;plot([-0.3 0.7],[0 0],'k')
        
        %p=0.0064 corr-coef non-directional
          %%% CSplot FOR SURROUNDING SYLLABLES
            % Open Adjsylls.mat and PredictwithZF.mat
          for i=1:12
              expnum=Adjsylls(i).expnum;
              if isequal(Predict(expnum).direction,'up')
                  prct=60;
              else
                  prct=40;
              end
              indices=ContingSimIND(Predict(expnum).Targeting,Adjsylls(i).pitchApre,prct);
              Prediction(i)=(mean(mean(Adjsylls(i).pitchBpre(Adjsylls(i).regB,indices)))-mean(mean(Adjsylls(i).pitchBpre(Adjsylls(i).regB,:))))/Adjsylls(i).Ashift;

          end
        figure;hold on;subplot(221);plot(Prediction,bac,'.','MarkerSize',20)
        corrcoef(Prediction,bac)
        subplot(222);plot(aac,bac,'.','MarkerSize',20)
        subplot(224);hold on;
        for i=1:12
            plot(Adjsylls(i).timeBafterA,Prediction(i),'.','MarkerSize',20,'Color','r')
            plot(Adjsylls(i).timeBafterA,bac(i),'.','MarkerSize',20)
        end
  % run 416.m  lines 1941 to 1950
  hold on;plot(t(201:902)/8,mnabb(358:1059)/j1,'b','LineWidth',2)
          plot(t(201:902)/8,mnabbT(358:1059)/j2,'r','LineWidth',2) % targ imprecision included


%%%%%% DIR Song
            Predict(15).ResidDIR=residspk30;
            Predict(18).ResidDIR=residspk28
            Predict(22).ResidDIR=residso11;
                        indDIR=[7 11 14 15 18 22]
            for i=indDIR
                    if isequal(Predict(i).direction,'up')
                    prctile=60;
                     else
                    prctile=40;
                    end
                    CSINA=mean(ContingSim2(Predict(i).Targeting,Predict(i).ResidINA,prctile));
                    if i<24 % long notes
                                        j=i
                                            count=0;
                                            for kk=0.5:0.01:2
                                                count=count+1;
                                                choices(count)=sum(abs((kk*abs(CSINA(Predict(i).onset:Predict(i).onset+350))./max(abs(CSINA(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                                            end
                                            disterINA(i)=min(choices);

                    else % short notes
                                        j=i
                                            count=0;
                                            for kk=0.5:0.01:2
                                                count=count+1;
                                                choices(count)=sum(abs((kk*abs(CSINA(Predict(i).onset:Predict(i).onset+160))./max(abs(CSINA(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                                            end
                                            disterINA(i)=min(choices);
                    end
            end
            indDIR=[7 11 14 15 18 22]
            for i=indDIR
                    if isequal(Predict(i).direction,'up')
                    prctile=60;
                     else
                    prctile=40;
                    end
                    CSDIR=mean(ContingSim2(Predict(i).Targeting,Predict(i).ResidDIR,prctile));
                    if i<24 % long notes
                                        j=i
                                            count=0;
                                            for kk=0.5:0.01:2
                                                count=count+1;
                                                choices(count)=sum(abs((kk*abs(CSDIR(Predict(i).onset:Predict(i).onset+350))./max(abs(CSDIR(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                                            end
                                            disterDIR(i)=min(choices);

                    else % short notes
                                        j=i
                                            count=0;
                                            for kk=0.5:0.01:2
                                                count=count+1;
                                                choices(count)=sum(abs((kk*abs(CDIR(Predict(i).onset:Predict(i).onset+160))./max(abs(CSDIR(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                                            end
                                            disterDIR(i)=min(choices);
                    end
            end
            for i=indDIR
                CVratio(i)=median(std(Predict(i).ResidDIR(Predict(i).onset:Predict(i).offset,:)')./std(Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:)'));
                %nAC(i)=round(CVratio(i)^2*size(Predictx(i).ResidINA,2));
                nAC(i)=size(Predict(i).ResidDIR,2);
                gamma=0.5;
                %nAC(i)=round(size(Predict(i).ResidDIR,2)*(CVratio(i)/(gamma+(1-gamma)*CVratio(i)))^2);
                if isequal(Predict(i).direction,'up')
                    prctile=60;
                else
                    prctile=40;
                end

                for k=1:size(Predict(i).ResidAC,2)-nAC(i)
                    sample=ceil(rand(nAC(i),1)*size(Predict(i).ResidAC,2));
                    CSguess=mean(ContingSim2(Predict(i).Targeting,Predict(i).ResidAC(:,sample),prctile));
                    CSguess=CSguess-mean(Predict(i).ResidAC(:,sample)');
                    CSguesses(k).data=CSguess;
                    % normalize and calculate distance
                    if i<24 % long notes
                                        j=i
                                            count=0;
                                            for kk=0.5:0.01:2
                                                count=count+1;
                                                choices(count)=sum(abs((kk*abs(CSguess(Predict(i).onset:Predict(i).onset+350))./max(abs(CSguess(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                                            end
                                            dister(i,k)=min(choices);

                    else % short notes
                                        j=i
                                            count=0;
                                            for kk=0.5:0.01:2
                                                count=count+1;
                                                choices(count)=sum(abs((kk*abs(CSguess(Predict(i).onset:Predict(i).onset+160))./max(abs(CSguess(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                                            end
                                            dister(i,k)=min(choices);
                    end
                end
            end
            median(dister')
            disterDIR
            % after this correction, AC does better than DIR for 2/3 cases and
            % slightly worse for 1/3

           % THIS SHOWS that like INA, DIR has more precise variation
            figure;
            for k=1:6
                i=indDIR(k)
                clear ina
                clear ac
                for j=1:size(Predict(i).ResidINA,2)
                    ina(j,:)=(xcorr(Predict(i).ResidINA(Predict(i).onset:Predict(i).onset+350,j)));
                    ina(j,:)=ina(j,:)/max(ina(j,:));
                end
                for j=1:size(Predict(i).ResidAC,2)
                    ac(j,:)=(xcorr(Predict(i).ResidAC(Predict(i).onset:Predict(i).onset+350,j)));
                    ac(j,:)=ac(j,:)/max(ac(j,:));
                end
                for j=1:size(Predict(i).ResidDIR,2)
                    dir(j,:)=(xcorr(Predict(i).ResidDIR(Predict(i).onset:Predict(i).onset+350,j)));
                    dir(j,:)=dir(j,:)/max(dir(j,:));
                end

                plot(mean(ac))
                hold on;plot(mean(ina),'r')
                hold on;plot(mean(dir),'k')
                %distbt(k)=mean(mean(ac)-mean(ina));
                % fwhm
                distbt(k)=(350-min(find(mean(ac)>0.5)))/(350-min(find(mean(ina)>0.5)));
            end
%%% CONSENSUS - BANDWIDTH etc.
   % ConsensusData structure
% Look at middle of Resid*4+512 and then take 600pts on either side!
        % Predict(25).data 
        % Alldata1sig(10) - exp 1(intact) vs. exp2 (lesion)
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=512+4*Alldata1sig(10).startnote-50;
                    last=512+4*Alldata1sig(10).endnote+50;
                    Intact(i,:)=ConsentPlot2(Alldata1sig(10).exp(1).rawdata(i,first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=512+4*Alldata1sig(10).startnote-50;
                    last=512+4*Alldata1sig(10).endnote+50;
                    Lesion(i,:)=ConsentPlot2(Alldata1sig(10).exp(2).rawdata(i+30,first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(7).data
        % Alldata1sig(12) - exp 1 (intact) vs. exp 2 stationary mean part (lesion)
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=2800;%512+4*Alldata1sig(12).startnote;
                    last=4000;%512+4*Alldata1sig(12).endnote;
                    Intact(i,:)=ConsentPlot2(Alldata1sig(12).exp(1).rawdata(i,first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=2800;%512+4*Alldata1sig(12).startnote;
                    last=4000%512+4*Alldata1sig(12).endnote;
                    Lesion(i,:)=ConsentPlot2(Alldata1sig(12).exp(2).rawdata(i+130,first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(26).data
        % Alldata1sig(11) - exp 1 (intact) vs. exp 2 (lesion)
                    Intact=[];
            Lesion=[];
                for i=1:30
                    first=512+4*Alldata1sig(11).startnote-100;
                    last=512+4*Alldata1sig(11).endnote+100;
                    Intact(i,:)=ConsentPlot2(Alldata1sig(11).exp(1).rawdata(i,first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=512+4*Alldata1sig(11).startnote-100;
                    last=512+4*Alldata1sig(11).endnote+100;
                    Lesion(i,:)=ConsentPlot2(Alldata1sig(11).exp(2).rawdata(i+10,first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(14)
        % bk63
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=2200;
                    last=3400;
                    Intact(i,:)=ConsentPlot2(fvals63AC1201(i).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=2200;
                    last=3400;
                    Lesion(i,:)=ConsentPlot2(fvals63MU1203(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(22)
        % o11bk33
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=3700;
                    last=4900;
                    Intact(i,:)=ConsentPlot2(fvalsPre(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=3700;
                    last=4900;
                    Lesion(i,:)=ConsentPlot2(fvalsPost(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(19)
        % bk13w63
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=1700;
                    last=2900;
                    Intact(i,:)=ConsentPlot2(fvalsPre(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=1700;
                    last=2900;
                    Lesion(i,:)=ConsentPlot2(fvalsPost(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(21)
        % o50bk72
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=2700;
                    last=3900;
                    Intact(i,:)=ConsentPlot2(fvalsPre(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=2700;
                    last=3900;
                    Lesion(i,:)=ConsentPlot2(fvalsPost(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(18)
        % pk28bk66
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=3500;
                    last=4700;
                    Intact(i,:)=ConsentPlot2(fvalsPre(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=3500;
                    last=4700;
                    Lesion(i,:)=ConsentPlot2(fvalsPost(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(16)
        % pk42bk73
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=2900;
                    last=4100;
                    Intact(i,:)=ConsentPlot2(fvalsPre(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=2900;
                    last=4100;
                    Lesion(i,:)=ConsentPlot2(fvalsPost(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(15)
        % pk30bk79
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=3600;
                    last=4800;
                    Intact(i,:)=ConsentPlot2(fvals1208(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=3600;
                    last=4800;
                    Lesion(i,:)=ConsentPlot2(fvals107(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(27)
        % pu57w52
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=900;
                    last=2000;
                    Intact(i,:)=ConsentPlot2(DataA(11).fvals(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=900;
                    last=2000;
                    Lesion(i,:)=ConsentPlBK63ot2(DataA(9).fvals(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
        % Predict(28)
        % pu56w26
            Intact=[];
            Lesion=[];
                for i=1:30
                    first=800;
                    last=1900;
                    Intact(i,:)=ConsentPlot2(fvals501ac(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
                for i=1:30
                    first=800;
                    last=1900;
                    Lesion(i,:)=ConsentPlot2(fvals430apv(i+10).datt(first:last),32000,1024,1020,1,3,5,5);
                end
medIntact=median([ConsensusData.P7intact;ConsensusData.P25intact;ConsensusData.P26intact;ConsensusData.P27intact;ConsensusData.P28intact;ConsensusData.P15intact;ConsensusData.P16intact;ConsensusData.P18intact;ConsensusData.P19intact;ConsensusData.P21intact;ConsensusData.P22intact]);
medLesion=median([ConsensusData.P7lesion;ConsensusData.P25lesion;ConsensusData.P26lesion;ConsensusData.P27lesion;ConsensusData.P28lesion;ConsensusData.P15lesion;ConsensusData.P16lesion;ConsensusData.P18lesion;ConsensusData.P19lesion;ConsensusData.P21lesion;ConsensusData.P22lesion]);

count=0;
vals=[0.1 0.3 0.5 0.7 0.9 1.1 1.3];
sInt=median(ConsensusData1.P22intact);
sLes=median(ConsensusData1.P22lesion);
[ymInt,xmInt]=max(sInt);
[ymLes,xmLes]=max(sLes);
count=count+1;
BandwidthInt(count)=pinterp([vals(xmInt-1);vals(xmInt);vals(xmInt+1)],[sInt(xmInt-1);sInt(xmInt);sInt(xmInt+1)]);
BandwidthLes(count)=pinterp([vals(xmLes-1);vals(xmLes);vals(xmLes+1)],[sLes(xmLes-1);sLes(xmLes);sLes(xmLes+1)]);
BandwidthDir(count)=pinterp([vals(xmLes-1);vals(xmLes);vals(xmLes+1)],[sLes(xmLes-1);sLes(xmLes);sLes(xmLes+1)]);
% notably BandwidthLes is a tiny bit smaller p=0.4230 - both are around 0.5

