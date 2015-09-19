function [fmrecentescapes,fmrecenthits,fmrecent]=trialbytrialcore(Experiment,expnumber,histwindow,whichwindow,comparisonmethod)
% run 'trialbytrialpre.m' beforehand
% []=trialbytrialcore(Experiment,100);
pitchWNon=Experiment(expnumber).pitchWNon;
pitchBaseline=Experiment(expnumber).pitchBaseline;
indEscapeAbove=Experiment(expnumber).indEscapeAbove;
indHitAbove=Experiment(expnumber).indHitAbove;
ptwindow=Experiment(expnumber).ptwindow;
if isequal(whichwindow,'targwindow');
    ptwindow=Experiment(expnumber).targwindow;
end
% PRE-PROCESSING
    % How much does pitch at t depend on pitch of the most recent n
    % performances. (n=histwindow)
    residWNon=[];
    for i=1:size(pitchWNon,2)
        residWNon(:,i)=pitchWNon(:,i)-mean(pitchBaseline')';
    end
    residWNon=residWNon(ptwindow,:); % only look at reasonable part

% REGRESSION
 b=zeros(size(residWNon,2)-histwindow,histwindow);
 
 
 if isequal(comparisonmethod,'subtraction')
        for i=histwindow+30:size(residWNon,2)-30
            mnpresent=mean(residWNon(:,i-29:i+30)');
            b(i-histwindow,1)=mean(abs(mnpresent-(residWNon(:,i-1)')));
            for j=2:histwindow
                x=mean(abs(mnpresent-mean(residWNon(:,i-j:i-1)')));
                b(i-histwindow,j)=x;
            end
        end  
        regcoefs=b;
 else
      if isequal(comparisonmethod,'division')
        for i=histwindow+1:size(residWNon,2)
            for j=1:histwindow
                x=median((residWNon(:,i)./residWNon(:,i-j)));
                b(i-histwindow,j)=x;
            end
        end  
        regcoefs=b;
      else
          if isequal(comparisonmethod,'corrcoef')
        for i=histwindow+1:size(residWNon,2)
                x=corrcoef(residWNon(:,i),(residWNon(:,i-1)));
                b(i-histwindow,1)=x(2);
            for j=2:histwindow
                x=corrcoef(residWNon(:,i),median(residWNon(:,i-j:i-1)'));
                b(i-histwindow,j)=x(2);  
            end
        end  
        regcoefs=b;
          end
      end
 end
 
% POST-PROCESSING
    for i=1:size(regcoefs,1)
        recentcoefs(i,:)=regcoefs(i,:);
        recentescapesabove(i,:)=indEscapeAbove(i:i+histwindow-1);
        recenthitsabove(i,:)=indHitAbove(i:i+histwindow-1);
    end


    for i=1:histwindow
        aind=find(recentescapesabove(:,i));
        if ~isempty(aind)
            holder=(recentcoefs(aind,i));
            mrecentescapes(i)=mean(holder); % Only look at the regression coefficients that are not outliers
        end
        bind=find(recenthitsabove(:,i));
        if ~isempty(bind)
            holder=(recentcoefs(bind,i));
            mrecenthits(i)=mean(holder); % Only look at the regression coefficients that are not outliers
        end
    %     abind=find(recentescapesabove(:,i)+recenthitsabove(:,i));
    %     if ~isempty(aind) | ~isempty(bind)
    %         mrecentabove(i)=median(recentcoefs(abind,i));
    %     end
    end
    mrecent=mean(recentcoefs);
  % Flip so that first point is most recent syllable.
  fmrecentescapes=fliplr(mrecentescapes);
  fmrecenthits=fliplr(mrecenthits);
  fmrecent=fliplr(mrecent);

