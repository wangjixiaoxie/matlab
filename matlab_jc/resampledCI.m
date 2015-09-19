function CIs=resampledCI(inputvec,alpha,ismean)
% alpha as 0.05
% ismean=1 - take mean;;; ismean=0 - take median
if ismean
    alphapercent=alpha*100;
    b=inputvec;
    for i=1:10000
        randassign=ceil(rand(1,length(inputvec))*length(inputvec));
        mm(i)=mean(b(randassign));
    end
    CIlow=prctile(mm,alphapercent);
    CIhigh=prctile(mm,100-alphapercent);
    CIs=[CIlow CIhigh];
else
        alphapercent=alpha*100;
    b=inputvec;
    for i=1:10000
        randassign=ceil(rand(1,length(inputvec))*length(inputvec));
        mm(i)=median(b(randassign));
    end
    CIlow=prctile(mm,alphapercent);
    CIhigh=prctile(mm,100-alphapercent);
    CIs=[CIlow CIhigh];
end
