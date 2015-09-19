function [phits,pescs,phitsNORM,pescsNORM,pescsHINORM,hits,escs]=TSA1014_2(vals)
fv=vals(:,2);
songs=vals(:,4);
meanmem=80;

for iii=1:20
    mem=iii;
posthitcount=0;
postesccount=0;
postescHIcount=0;
for i=mem+1:length(fv)

        if vals(i-mem,6)==1 %If pvs note was a wnblast
            posthitcount=posthitcount+1;
            hits(posthitcount)=fv(i-mem);
            phits(posthitcount)=fv(i);
            songA=songs(i);
            nextsong=[];
            count=0;
            k=i-mem;
            while 1
                k=k+1;
                if k==length(songs) break; end
                if songs(k)~=songA
                    if count>0 && songs(k)~=songB
                        break
                    end
                    count=count+1;
                    nextsong(count)=fv(k);
                    songB=songs(k);
                end
            end
            if nextsong>2000
            phitsNEXT(posthitcount)=mean(nextsong);
            end
        else %If pvs note was an escape
            postesccount=postesccount+1;
            escs(postesccount)=fv(i-mem);
            pescs(postesccount)=fv(i);
            songA=songs(i);
            nextsong=[];
            count=0;
            k=i-mem;
            while 1
                k=k+1;
                if k==length(songs) break; end
                if songs(k)~=songA
                    if count>0 && songs(k)~=songB
                        break
                    end
                    count=count+1;
                    nextsong(count)=fv(k);
                    songB=songs(k);
                end
            end
            if nextsong>2000
            pescsNEXT(postesccount)=mean(nextsong);
            end
        end
end
        meanpescs(iii)=mean(pescs);
        meanphits(iii)=mean(phits);
end
posthitcount=0;
postesccount=0;
postescHIcount=0;
for i=meanmem+1:length(fv)
    last50mean=mean(fv(i-meanmem:i-mem));
    fvNORM=fv(i)-last50mean;
    if vals(i-mem,5)==1 %If pvs note was below threshold
        if vals(i-mem,6)==1 %If pvs note was a wnblast
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
        