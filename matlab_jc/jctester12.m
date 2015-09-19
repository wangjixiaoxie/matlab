function B=jctester12(ACSF,ontime,offtime) 
for i=1:length(ACSF)
        FFF=[];
        count=0;
        clear pitch;
        
            pitch=ACSF(i).pitchHILBERT4;
        for k=1:size(pitch,2)
            dept=0;
            a=pitch(ontime:offtime,k);
            for j=11:length(a)
                b(j-10)=abs(a(j)-a(j-10));
            end
            if max(b)>30
                dept=1;
            end
            if dept==0
                count=count+1;
                FFF(:,count)=pitch(:,k);
            end
        end
        B(i).selectedpitchcurves=FFF;
end

