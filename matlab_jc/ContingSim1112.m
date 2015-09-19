function selectedcurves=ContingSim1112(targpoint,Alldata,cont)

for k=3:length(Alldata)
    pitch=Alldata(k).baselineAC(Alldata(k).startnote:Alldata(k).endnote,:);
    mnbase=mean(pitch');
    for j=1:size(pitch,2) % for each note
        pointvalues(j)=pitch(targpoint,j);
    end
        contingency=prctile(pointvalues,cont);
        count=1;
        clear selectedcurves
    for i=1:size(pitch,2)
        if pitch(targpoint,i)>=contingency
            selectedcurves(:,count)=pitch(:,i);
            count=count+1;
        end
    end
    hold on;plot(mean(selectedcurves')-mnbase)
end

%for k=1:20
%     for i=1:size(pitch,2)
%         rr=round(rand*length(toffsets));
%         if rr==0
%             rr=1;
%         end
%         offset=round(toffsets(1)); %round(mean(toffsets)+0.8*std(toffsets)*randn);
%         ffs(i)=mean(pitch(offset-32:offset,i)); % estimate of pitch within the window
%     end
%     L=prctile(ffs,cont);
%     ind=find(ffs<L);
%     aav=zeros(size(pitch,1),1);
%     for j=1:length(ind)
%         aav=aav+pitch(:,ind(j));
%     end
%     aavTot=aav/j;
%end
%aavfin=mean(aavTot);