function [meanphits,meanpescs,medianphits,medianpescs,meanNextesc,meanNexthit,hits,escs]=TSA1014_3(vals)
%same position in next song
fv=vals(:,2);
songs=vals(:,4);
meanmem=80;

NotePositions(1)=1;
Nowsong=songs(1);
for ii=2:length(songs)
    if songs(ii)==Nowsong
        NotePositions(ii)=NotePositions(ii-1)+1;
    else
        NotePositions(ii)=1;
        Nowsong=songs(ii);
    end
end




for iii=1:20
    mem=iii;
posthitcount=0;
postesccount=0;

for i=mem+1:length(fv)-20
        NotePosit=NotePositions(i);
        gg=NotePositions(i+1:length(NotePositions));
        hh=find(gg==NotePosit);
        if length(hh)<3 || hh(3)+i+1>length(fv)
            NextSongSpot=0;
        else
        NextSongSpot=fv(hh(3)+i);  %%%%adding 1 gives no relation!!!!
        end
        if vals(i-mem,6)==1 %If pvs note was a wnblast
            posthitcount=posthitcount+1;
            hits(posthitcount)=fv(i-mem);
            phits(posthitcount)=fv(i);
            
            Nextsonghit(posthitcount)=NextSongSpot;
        else %If pvs note was an escape
            postesccount=postesccount+1;
            escs(postesccount)=fv(i-mem);
            pescs(postesccount)=fv(i);
            Nextsongesc(postesccount)=NextSongSpot;
        end
end
        if iii==1
            meanNextesc=mean(Nextsongesc(find(Nextsongesc>1)));
            meanNexthit=mean(Nextsonghit(find(Nextsonghit>1)));
        end
        meanpescs(iii)=mean(pescs);
        medianpescs(iii)=median(pescs);
        meanphits(iii)=mean(phits);
        medianphits(iii)=median(phits);
end
% posthitcount=0;
% postesccount=0;
% postescHIcount=0;
% for i=meanmem+1:length(fv)
%     last50mean=mean(fv(i-meanmem:i-mem));
%     fvNORM=fv(i)-last50mean;
%     if vals(i-mem,5)==1 %If pvs note was below threshold
%         if vals(i-mem,6)==1 %If pvs note was a wnblast
%             posthitcount=posthitcount+1;
%             phitsNORM(posthitcount)=fvNORM;
%         else
%             postesccount=postesccount+1;
%             pescsNORM(postesccount)=fvNORM;
%         end
%     else
%         postescHIcount=postescHIcount+1;
%         pescsHINORM(postescHIcount)=fvNORM;
%     end
% end
%         