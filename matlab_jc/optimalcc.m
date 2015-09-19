function [bestcrossco,mindist]=optimalcc(exp,AC)    
% [Predict(1).bestcrosscoAC,Predict(1).mindistAC]=optimalcc(Predict(1),1);
% 
if AC==1
    default=1./max(exp.crosscoAC(exp.onset:exp.offset));
        for i=1:200
            fac(i)=default-0.8+0.005*i;
            matchscore(i)=sum(abs(exp.crosscoAC(exp.onset:exp.offset)*fac(i)-exp.LearnedNorm(exp.onset:exp.offset)));
        end
    [mindist,ind]=min(matchscore);
    bestcrossco=exp.crosscoAC*fac(ind);
else
    default=1./max(exp.crosscoINA(exp.onset:exp.offset));
        for i=1:200
            fac(i)=default-0.8+0.005*i;
            matchscore(i)=sum(abs(exp.crosscoINA(exp.onset:exp.offset)*fac(i)-exp.LearnedNorm(exp.onset:exp.offset)));
        end
    [mindist,ind]=min(matchscore);
    bestcrossco=exp.crosscoINA*fac(ind);
end

