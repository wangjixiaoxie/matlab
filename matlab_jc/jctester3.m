function B=jctester3(B,ontime,offtime) 
for i=1:length(B)
        FFF=[];
        count=0;
        clear pitch;
        
            pitch=B(i).rawpitchcurves;
        for k=1:size(pitch,2)
            dept=0;
            a=pitch(ontime:offtime,k);
            for j=2:length(a)
                b(j-1)=abs(a(j)-a(j-1));
            end
            if max(b)>50
                dept=1;
            end
            if dept==0
                count=count+1;
                FFF(:,count)=pitch(:,k);
            end
        end
        B(i).selectedpitchcurves=FFF;
end

