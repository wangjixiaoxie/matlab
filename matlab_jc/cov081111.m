load /bulbul2/CovertAnalysis/Covert040810.mat
Experiment=[Exps1 Exps2 Exps3];
figure;hold on;
thistime=[250 250 250 200 200 200 300 300 300 300 300 250 350 350 250 300 300 250 250 300];

clear mnSEpost mnSEpre1 mnSEpre2 mnAMPpre1 mnAMPpre2 mnAMPpost 
indX=ind([1:14 16:end]);
for k=1:length(indX)
    i=indX(k)
    timechunk=thistime(k)*4+512-128:thistime(k)*4+512+128;
    Fs=32000;
    NFFT = 2^nextpow2(length(timechunk)); % Next power of 2 from length of y
    sdpdfPRE=[];ampPRE=[];clear tPRE
    Epre=Experiment(i).fvACpre;
    for ii=1:length(Epre)
        ampPRE(ii)=sum(abs(Epre(ii).datt(timechunk)));
        Y = fft(Epre(ii).datt(timechunk)-mean(Epre(ii).datt(timechunk)),NFFT)/length(timechunk);
        sd=abs(Y(1:NFFT/2+1));
        sdpdfPRE(ii,:)=sd/sum(sd);
        durPRE(ii)=Epre(ii).offs(Epre(ii).ind)-Epre(ii).ons(Epre(ii).ind);
        tPRE(ii)=timing3(Epre(ii));
    end
    sdpdfPOST=[];ampPOST=[];clear tPOST
    Epost=Experiment(i).fvACpost(postday(i,1):postday(i,2));
    for ii=1:length(Epost)
        ampPOST(ii)=sum(abs(Epost(ii).datt(timechunk)));
        Y = fft(Epost(ii).datt(timechunk)-mean(Epost(ii).datt(timechunk)),NFFT)/length(timechunk);
        sd=abs(Y(1:NFFT/2+1));
        sdpdfPOST(ii,:)=sd/sum(sd);
        durPOST(ii)=Epost(ii).offs(Epost(ii).ind)-Epost(ii).ons(Epost(ii).ind);
        tPOST(ii)=timing3(Epost(ii));
    end
    %octavespacing=round(mean(E(thistime(k),:))/(Fs/2/256))/2;
    [a b]=max(mean(sdpdfPRE));
    clear spectentropyPRE
    for ii=1:size(sdpdfPRE,1)
        dat=sdpdfPRE(ii,:);%b-octavespacing:b+octavespacing);
        spectentropyPRE(ii)=-1*sum(dat.*log(dat));
    end
    clear spectentropyPOST
    for ii=1:size(sdpdfPOST,1)
        dat=sdpdfPOST(ii,:);%b-octavespacing:b+octavespacing);
        spectentropyPOST(ii)=-1*sum(dat.*log(dat));
    end
    
    % duration
    
   tdelt1(k)=(tPRE(round(length(tPRE)/2)))-(tPRE(1:))
   tdelt2(k)=mean(tPOST)-mean(tPRE(3*round(length(tPRE)/4):end));
   realsPRE=spectentropyPRE(~isnan(spectentropyPRE));
   mnSEpre1(k)=mean(realsPRE(1:round(length(realsPRE)/2)));
   mnSEpost(k)=mean(spectentropyPOST(~isnan(spectentropyPOST)));
   mnSEpre2(k)=mean(realsPRE(round(length(realsPRE)/2):end));
   mnAMPpre1(k)=mean(ampPRE(1:round(length(ampPRE)/2)));
   mnAMPpre2(k)=mean(ampPRE(round(length(ampPRE)/2):end));
   mnAMPpost(k)=mean(ampPOST);
   mnDURpre1(k)=mean(durPRE(1:round(length(durPRE)/2)));
   mnDURpre2(k)=mean(durPRE(round(length(durPRE)/2):end));
   mnDURpost(k)=mean(durPOST);
   
   
   sdSEpre1(k)=std(realsPRE(1:round(length(realsPRE)/2)));
   sdSEpost(k)=std(spectentropyPOST(~isnan(spectentropyPOST)));
   sdSEpre2(k)=std(realsPRE(round(length(realsPRE)/2):end));
   sdAMPpre1(k)=std(ampPRE(1:round(length(ampPRE)/2)));
   sdAMPpre2(k)=std(ampPRE(round(length(ampPRE)/2):end));
   sdAMPpost(k)=std(ampPOST);
   sdDURpre1(k)=std(durPRE(1:round(length(durPRE)/2)));
   sdDURpre2(k)=std(durPRE(round(length(durPRE)/2):end));
   sdDURpost(k)=std(durPOST);
  
   
   mnFFpre1(k)=mean(mean(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,1:round(size(Experiment(i).pitchACpre,2)/2))));
   mnFFpre2(k)=mean(mean(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,round(size(Experiment(i).pitchACpre,2)/2):end)));
   mnFFpost(k)=mean(mean(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,:)));
   sdFFpre1(k)=std(mean(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,1:round(size(Experiment(i).pitchACpre,2)/2))));
   sdFFpre2(k)=std(mean(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,round(size(Experiment(i).pitchACpre,2)/2):end)));
   sdFFpost(k)=std(mean(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,:)));
   
    coef(k)=1-2*(isequal(Experiment(i).DIR,'down'));
end

figure;hold on;
subplot(121);hold on;
plot(mnSEpre1,mnSEpre2,'.');
plot(mnSEpre2,mnSEpost,'r.');
xlim([4 5]);ylim([4 5]);plot([4 5],[4 5]) % p=0.78

subplot(122);hold on;
plot(log(mnAMPpre1),log(mnAMPpre2),'.');
plot(log(mnAMPpre2),log(mnAMPpost),'r.'); % p=0.18
xlim([9 15]);ylim([9 15]);plot([9 15],[9 15]) 

figure;hold on;
subplot(121);hold on;
plot(mnSEpre1-mnSEpre2,mnSEpost-mnSEpre2,'k.')
plot([-0.2 0.1],[-0.2 0.1]);xlim([-0.2 0.1]);ylim([-0.2 0.1]) % p=0.65 sign test
subplot(122);hold on;
plot(log(mnAMPpre1([1:7 9:end]))-log(mnAMPpre2([1:7 9:end])),log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end])),'k.') % p=0.24 signtest
plot([-0.5 1],[-0.5 1]);xlim([-0.5 1]);ylim([-0.5 1])

subplot(133);hold on;
plot(coef.*(mnFFpre1-mnFFpre2),coef.*(mnFFpost-mnFFpre2),'k.') % p=0.02
plot([-110 50],[-110 50]);xlim([-110 50]);ylim([-110 50])

figure;hold on;
bar(1,mean(mnSEpre1-mnSEpre2))
errorbar(1,mean(mnSEpre1-mnSEpre2),std(mnSEpre1-mnSEpre2)/sqrt(length(mnSEpre1)))
bar(2,mean(mnSEpost-mnSEpre2))
errorbar(2,mean(mnSEpost-mnSEpre2),std(mnSEpost-mnSEpre2)/sqrt(length(mnSEpost)))
bar(3,mean(log(mnAMPpre1([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))
errorbar(3,mean(log(mnAMPpre1([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))),std(log(mnAMPpre1([1:7 9:end]))-log(mnAMPpre2([1:7 9:end])))/sqrt(length(log(mnAMPpre1([1:7 9:end])))))
bar(4,mean(log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))
errorbar(4,mean(log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))),std(log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end])))/sqrt(length(log(mnAMPpost([1:7 9:end])))))

% This is best
    figure;hold on;
    bar(1,mean(abs(mnSEpre1-mnSEpre2))/mean(sdSEpre2))
    errorbar(1,mean(abs(mnSEpre1-mnSEpre2))/mean(sdSEpre2),std(abs(mnSEpre1-mnSEpre2))/sqrt(length(mnSEpre1))/mean(sdSEpre2))
    bar(2,mean(abs(mnSEpost-mnSEpre2))/mean(sdSEpre2))
    errorbar(2,mean(abs(mnSEpost-mnSEpre2))/mean(sdSEpre2),std(abs(mnSEpost-mnSEpre2))/sqrt(length(mnSEpost))/mean(sdSEpre2))

    bar(4,mean(abs(log((mnAMPpre1([1:7 9:end])))-log((mnAMPpre2([1:7 9:end])))))/mean(log(sdAMPpre2([1:7 9:end]))))
    errorbar(4,mean(abs(log(mnAMPpre1([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))/mean(log(sdAMPpre2)),std(abs(log(mnAMPpre1([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))/sqrt(length(log(mnAMPpre1([1:7 9:end]))))/mean(log(sdAMPpre2([1:7 9:end]))))
    bar(5,mean(abs(log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))/mean(log(sdAMPpre2([1:7 9:end]))))
    errorbar(5,mean(abs(log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))/mean(log(sdAMPpre2([1:7 9:end]))),std(abs(log(mnAMPpost([1:7 9:end]))-log(mnAMPpre2([1:7 9:end]))))/sqrt(length(log(mnAMPpost([1:7 9:end]))))/mean(log(sdAMPpre2([1:7 9:end]))))

    bar(7,mean(abs(mnDURpre1-mnDURpre2))/mean(sdDURpre2))
    errorbar(7,mean(abs(mnDURpre1-mnDURpre2))/mean(sdDURpre2),std(abs(mnDURpre1-mnDURpre2))/sqrt(length(mnDURpre1))/mean(sdDURpre2))
    bar(8,mean(abs(mnDURpost-mnDURpre2))/mean(sdDURpre2))
    errorbar(8,mean(abs(mnDURpost-mnDURpre2))/mean(sdDURpre2),std(abs(mnDURpost-mnDURpre2))/sqrt(length(mnDURpost))/mean(sdDURpre2))

%   bar(7,mean(coef.*(mnFFpre1-mnFFpre2))/mean(sdFFpre2))
% errorbar(7,mean(coef.*(mnFFpre1-mnFFpre2))/mean(sdFFpre2),std(coef.*(mnFFpre1-mnFFpre2))/sqrt(length(mnFFpre1))/mean(sdFFpre2))
% bar(8,mean(Learning(indX))/mean(sdFFpre2))
% errorbar(8,mean((Learning(indX)))/mean(sdFFpre2),std((Learning(indX)))/sqrt(length(mnFFpost))/mean(sdFFpre2))  
% Do it for positive control experiments --- but I don't have raw data
        figure;hold on;
        thistime=[250 300 280 250 250 300 350 300 350 250 300 350 300 280];
