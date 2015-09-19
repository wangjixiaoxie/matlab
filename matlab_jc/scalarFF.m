function scalarFF(indE,medianpitch)
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=151:length(medianpitch)
    res=medianpitch(i)-median(medianpitch(i-15:i-10));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            %Hit=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=res;
            
        else
            %Esc=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=res;
        end
    end
end
for ii=1:size(facHit,2)
    g=min(abs(facHit(:,ii)));
    if g==0
        sizeHit=ii-1;
        break
    end
end
for ii=1:size(facEsc,2)
    g=min(abs(facEsc(:,ii)));
    if g==0
        sizeEsc=ii-1;
        break
    end
end
figure;plot(median(facHit(:,1:sizeHit)'))
hold on;plot(median(facEsc(:,1:sizeEsc)'),'r')
hold on;plot([0 50],[0 0],'k')
hold on;plot(median(facHit(:,1:sizeHit)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:sizeEsc)'),'.','MarkerSize',15,'Color','r')