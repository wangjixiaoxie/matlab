% TSA720
pitchTSAhar1=jc_pitchmat1024(shiftedTSA,1024,1020,1,1800,2700,[1],'obs0',1);
figure;plot(median(pitchTSAhar1(390:425,:)),'*')
% Look for a bimodal distribution
figure;hist(std(pitchTSAhar1(550:700,:)))
indE=find((std(pitchTSAhar1(550:700,:))<100));
figure;plot(pitchTSAhar1(:,indE))

%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:size(pitchTSAhar1,2)
    res=median((pitchTSAhar1(420:460,i)-median(pitchTSAhar1(420:460,i-20:i-1)')')');
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            Hit=median((pitchTSAhar1(420:460,index)-median(pitchTSAhar1(420:460,i-20:i-1)')')');
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=res-Hit;
            
        else
            Esc=median((pitchTSAhar1(420:460,index)-median(pitchTSAhar1(420:460,i-20:i-1)')')');
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=res-Esc;
        end
    end
end
figure;plot(median(facHit(:,1:199)'))
hold on;plot(median(facEsc(:,1:78)'),'r')
hold on;plot([0 100],[0 0],'k')
hold on;plot(median(facHit(:,1:199)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:78)'),'.','MarkerSize',15,'Color','r')