function [pretime posttime]=ContingSimMSB(toffsets,pitch,cont,thresholdCriterion)
%aavfin=ContingSimMSB(Predict(1).Targeting-Predict(1).onset,Predict(1).ResidAC(Predict(1).onset:Predict(1).onset+350,:),60);
% pitch - just the residuals within the note




% threshold criterion - A low FF threshold indicates that you are at the
% contingency position.  A high FF threshold indicates you are not near it.
% This value should be around 0.8 to provide the most information.  This
% can be set to whatever is optimal for all experiments.
%thresholdCriterion=0.8;

% decision criterion - escape rate that satisfies the bird - should be about the same as hit rate at asymptote
decisionCriterion=0.8;

% How many notes does the bird really get to use? Assume 8000 trials.
% Sometimes we have more info than this.  Is this a problem?


pitch=pitch(:,floor(rand(1,50)*size(pitch,2)+1));


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
        if cont>50
            ind=find(ffs>L);
        else
            ind=find(ffs<L);
        end
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
% for i=1:size(pitch,2)
%     for j=1:size(pitch,1)
%         FFvalue(j,i)=pitch(j,i);
%         EscapeLikelihood(j,i)=histog(j);
%     end
% end
% bin them\

% This loop does the transformation from curves that escape to points that
% escape
for i=1:size(pitch,1)
    FF=-0.1:0.001:0.1;
    for j=1:201
        ind1=find(pitch(i,:)>FF(j));
        if ~isempty(ind1)
            ind2=find(pitch(i,ind1)<FF(j)+0.01);
            if ~isempty(ind2)
                Binned(i,j)=mean(histog(ind1(ind2)));
            else
                Binned(i,j)=0;
            end
        else
            Binned(i,j)=0;
        end
    end
end
% get threshold crossings - max if 40, min 
for i=1:size(pitch,1)
    if max(Binned(i,1:201)>thresholdCriterion)
        thrcross(i)=FF(min(find((Binned(i,1:201)>thresholdCriterion))));
    else
        thrcross(i)=FF(201);
    end
end

% Decision: start at most likely 5ms and move outward until decision
% criterion is met
wideness=40; % 5ms
    for i=1:size(pitch,1)-wideness
        thr5ms(i)=median(thrcross(i:i+wideness-1));
    end
    [a,centershift]=min(thr5ms);
    centershift=centershift+wideness/2;
    % How far outward do you have to move?
    for i=1:size(pitch,1)
        if centershift-i<=0
            startpt(i)=1;
        else
            startpt(i)=centershift-i;
        end
        if centershift+i>=size(pitch,1)
            endpt(i)=size(pitch,1);
        else
            endpt(i)=centershift+i;
        end
        ind1=find(toffsets<endpt(i));
        ind2=find(toffsets(ind1)>startpt(i));
        targsaccounted(i)=length(ind2)/length(toffsets);
    end


    decisionCriterion=0.8;
    valu=min(find(targsaccounted>decisionCriterion));
    pretime=startpt(valu)-median(toffsets);
    posttime=endpt(valu)-median(toffsets);
    

g=7;
