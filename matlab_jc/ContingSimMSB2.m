function FFnorm=ContingSimMSB2(toffsets,pitch,cont,direction,minstdGauss,heightShift)
%aavfin=ContingSimMSB(Predict(1).Targeting-Predict(1).onset,Predict(1).ResidAC(Predict(1).onset:Predict(1).onset+350,:),60);
% pitch - just the residuals within the note


    % FFnorms(i,:)=ContingSimMSB2(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),80,Predict(i).direction,30,80);

% threshold criterion - A low FF threshold indicates that you are at the
% contingency position.  A high FF threshold indicates you are not near it.
% This value should be around 0.8 to provide the most information.  This
% can be set to whatever is optimal for all experiments.
%thresholdCriterion=0.8;

% decision criterion - escape rate that satisfies the bird - should be about the same as hit rate at asymptote
decisionCriterion=0.8;


% A miserly approach is to change FF only where FF is a good predictor of
    % outcome (outcome meaning WN vs. no WN).  Thus we ask: where does FF
    % predict outcome the best?
    
    % To quantify outcome, we simulate to determine the probability that
    % each curve will escape.  The bird could accomplish this by performing 
    % many more notes than we have baseline note and then making a
    % histogram.  For each time point, we do vector multiplication of the
    % FF vector (across curves) and the p(escape) vector (across curves).
   
        % How many notes does the bird really get to use? Assume 8000 trials.
        % Sometimes we have more info than this.  Is this a problem?
        pitch=pitch(:,floor(rand(1,50)*size(pitch,2)+1));
        
      % NORMALIZE TO UPSHIFT
        if ~isequal(direction,'up')
            pitch=pitch*-1;
        end
        
        % only consider targeting within the syllable
        ind1=find(toffsets>1);
        ind2=find(toffsets(ind1)<size(pitch,1));
        toffsets=toffsets(ind1(ind2));
        indices=[];
        count=0;
            for i=1:length(toffsets)
                offset=round(toffsets(i)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
                ffs=(pitch(offset,:)); % estimate of pitch within the window
                L=prctile(ffs,cont);
                ind=find(ffs>L);
                if ~isempty(ind)
                    count=count+1;
                    indices=[indices ind];
                end
            end
        aavfin=indices;
        % What proportion of the time does each curve escape?
        for i=1:size(pitch,2)
            histog(i)=length(find(aavfin==i))/length(toffsets);
        end
   %%%%%%%
   %%%%%%%
   %%%%%%%
   % This determines the 5ms (40pt) time region that is most predictive. 
    Mostpredictive=pitch*histog';
    for k=21:length(Mostpredictive)-20
        Mp5ms(k)=mean(Mostpredictive(k-20:k+20));
    end
    [a,indmaxMN]=max(Mp5ms);
% Now we build a FF Gaussian at the most predictive time region
    %minstdGauss=10; % width of motor primitive in ms
    heightGauss=prctile(pitch(indmaxMN,:),heightShift); % maximal shift height
    meanGauss=indmaxMN;
    FFshift=heightGauss*gaussian([1:size(pitch,1)],meanGauss,minstdGauss*8);
   % Find the time points in the Gaussian that clear threshold
    Gaussclears=find(FFshift>prctile(pitch(indmaxMN,:),cont));
   % What proportion of all targeting points included
        if isempty(Gaussclears)
            propcovered=0;
        else
            propcovered=sum(toffsets>min(Gaussclears-1) & toffsets < max(Gaussclears+1))/length(toffsets);
        end
% Iterative process - if this doesn't work, then we expand.  How do we
% expand?  Center at the best place.
    centerwidth=0;
    while propcovered < decisionCriterion
        centerwidth=centerwidth+16; % expand by 2ms
       % Where is the best place to have the center width?
        Mp5ms=[];
        for k=centerwidth/2+1:length(Mostpredictive)-centerwidth/2
            Mp5ms(k)=mean(Mostpredictive(k-centerwidth/2:k+centerwidth/2));
        end
        [a,indmax]=max(Mp5ms);
       % Put center width there and Gaussian decay on either side
        %minstdGauss=10; % width of motor primitive in ms
        heightGauss=prctile(pitch(indmax,:),heightShift); % maximal shift height
        meanGauss=indmax-centerwidth/2;
        if centerwidth>size(pitch,1)-3
            break
        end
        Gaussest=heightGauss*gaussian([1:size(pitch,1)-centerwidth],meanGauss,minstdGauss*8);
        FFshift=[Gaussest(1:meanGauss) heightGauss*(ones(1,centerwidth)) Gaussest(meanGauss+1:end)];
        Gaussclears=find(FFshift>prctile(pitch(indmaxMN,:),cont));
        if isempty(Gaussclears)
            propcovered=0;
        else
            propcovered=sum(toffsets>min(Gaussclears-1) & toffsets < max(Gaussclears+1))/length(toffsets);
        end
    end
    
    FFshift=FFshift/max(FFshift);
FFnorm=-1*ones(1,1500);
FFnorm(501-round(median(toffsets)):501-round(median(toffsets))+length(FFshift)-1)=FFshift;
