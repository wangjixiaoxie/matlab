function pitchavg=plotharmonicsson(data,N,OVERLAP,sigma,F_low,F_high,harmonics)
figure; hold on;

for i=1:harmonics
    pitches=jc_PitchData610son(data,N,OVERLAP,sigma,F_low,F_high,i);
    if i==1; str='red'; 
    else
        if i==2; 
            str='blue'; 
        else
            if i==3; 
                str='black'; 
            else
                if i==4; 
                    str='yellow';
                end
            end
        end
    end
    plot(pitches,'Color',str)
end
pitchavg=jc_PitchData610sonavg(data,N,OVERLAP,sigma,F_low,F_high,harmonics);
plot(pitchavg,'Color','green')
h = legend('harmonic1','harmonic2','harmonic3','avg',4);