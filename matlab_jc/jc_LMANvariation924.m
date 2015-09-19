function LMANcurves=jc_LMANvariation3(x1,x2)
% e.g. x1=ACSF(3).residHILB2(1900:2500,21:50);
% e.g. x2=INA(3).residHILB2(2000:2600,1:30);

%%%%%%%%
% The idea of this is two-fold.
% 1. produce the autocorrelation of the "LMAN-dependent variability"
% 2. simulate LMAN variability -- decompose signals containing and lacking
% LMAN, subtract the spectral difference and add that back to the signal
% containing LMAN --- TO BE USED for learning simulation




%%%%% 12-15-08 %%%%%%%%
    % ********normalize after taking the mean so that you take into account the
    % importance of the notes with big residuals
figure;plot(mean(jc_xcorr(x1)')./max(mean(jc_xcorr(x1)')),'r')
hold on;plot(mean(jc_xcorr(x2)')./max(mean(jc_xcorr(x2)')),'b')
plotLMANxcorr(x1,x2);

% This code plots the widest it can be
for i=1:size(x1,1)
    ll(i)=10;
end
hold on;plot(xcorr(ll)./max(xcorr(ll)),'k')

%%% Basic method 
    % This takes noLMAN (x2) and adds randomized magnitude differences between 
    % yesLMAN (x1) and noLMAN  at randomized phases selected from the
    % yesLMAN phase distribution (why this distn?).  This sum should approximately mimic the
    % timescale of yesLMAN variation, as confirmed by similarity in the
    % autocorrelation function.
clear LMANcurves
clear varLMAN
clear mag1
clear mag2 
clear pha1
clear pha2
long=min([size(x1,2) size(x2,2)]);
for i=1:long
    a=x1(:,i); 
    b=x2(:,i); 
    z1=fft(a);
    z2=fft(b);
    mag1(i,:)=abs(z1);
    mag2(i,:)=abs(z2);
    pha1(i,:)=angle(z1);
    pha2(i,:)=angle(z2);
end
% meanLMAN=mean(mag1)-mean(mag2);
% % variability from the median
% for i=1:size(mag1,1)
%     varLMAN(i,:)=(mag1(i,:)./mean(mag1)).*meanLMAN;
% end
% 
% %Simulation;
clear LMANcurves
count=0;

    for i=1:100
        count=count+1;
        m=median(mag2-mag1)+randn*std(mag2-mag1);
        %m=varLMAN(round(rand*(long-1))+1,:);
       
        %p=(pha1(round(rand*(long-1))+1,:)); % this could be changed.
        p=median(pha1)+randn*std(pha1)/5;
        re=m.*cos(p);     
        im=m.*sin(p);
        dat=complex(re,im);
        gg=real(ifft(dat));
        LMANcurves(:,count)=gg;
    end


            notelength=sizer;
            numms=floor(notelength/8);
            clear crosscoAC
            clear crosscoINA
            clear mnccAC
            clear mnccINA
            clear mnccLMAN
            clear crosscoLMAN
            for ii=1:numms
                first=ii*8;
                middle=firstpoint+first;
                init=500-first+1;
                for j=1:notelength
                    ab=corrcoef(x1(j,:),x1(middle-firstpoint,:));
                    crosscoAC(init+j,ii)=ab(2);
                    ai=corrcoef(x2(j,:),x2(middle-firstpoint,:));
                    crosscoINA(init+j,ii)=ai(2);
                    index1=find(~isnan(mean(signal)));
                    aL=corrcoef(LMANcurves(j,index1),LMANcurves(middle-firstpoint,index1));
                    crosscoLMAN(init+j,ii)=aL(2);
                end
            end
            for i=1:size(crosscoAC,1)
                ind1=find(crosscoAC(i,:)>0);
                mnccAC(i)=mean(crosscoAC(i,ind1));
                ind2=find(crosscoINA(i,:)>0);
                mnccINA(i)=mean(crosscoINA(i,ind2));
                ind3=find(crosscoLMAN(i,:)>0);
                mnccLMAN(i)=mean(crosscoLMAN(i,ind3));
            end

            figure;plot(mnccAC)
            hold on;plot(mnccINA,'r')
            plot(mnccLMAN,'g')
        


