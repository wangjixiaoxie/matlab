function [phits,pescs,pescsHI,phitsNORM,pescsNORM,pescsHINORM,hits,escs]=TSA1013(vals)
fv=vals(:,2);
posthitcount=0;
postesccount=0;
postescHIcount=0;
for i=5:length(fv)
    if vals(i-1,5)==1 %If pvs note was below threshold
        if vals(i-1,6)==1 %If pvs note was a wnblast
            posthitcount=posthitcount+1;
            hits(posthitcount)=fv(i-1);
            phits(posthitcount)=fv(i);
        else
            postesccount=postesccount+1;
            escs(postesccount)=fv(i-1);
            pescs(postesccount)=fv(i);
        end
    else
        postescHIcount=postescHIcount+1;
        pescsHI(postescHIcount)=fv(i);
    end
end

posthitcount=0;
postesccount=0;
postescHIcount=0;
for i=51:length(fv)
    last50mean=mean(fv(i-50:i-1));
    fvNORM=fv(i)-last50mean;
    if vals(i-1,5)==1 %If pvs note was below threshold
        if vals(i-1,6)==1 %If pvs note was a wnblast
            posthitcount=posthitcount+1;
            phitsNORM(posthitcount)=fvNORM;
        else
            postesccount=postesccount+1;
            pescsNORM(postesccount)=fvNORM;
        end
    else
        postescHIcount=postescHIcount+1;
        pescsHINORM(postescHIcount)=fvNORM;
    end
end
        