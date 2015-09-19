function B=jctester3(B,ontime,offtime) 
for i=1:length(B)
        FFF=[];
        count=0;
        %clear pitch;
        clear b;
          %pitch=Alldata(4).pitchDpost;
        for k=1:50
            dept=0;
            a=pitch(700:1200,k);
            for j=2:length(a)
                b(j-1)=abs(a(j)-a(j-1));
            end
            if max(b)>70
                dept=1;
            end
            if dept==0
                count=count+1;
                FFF(:,count)=pitch(:,k);
            end
        end
        B(i).selectedpitchcurves=FFF;
end

