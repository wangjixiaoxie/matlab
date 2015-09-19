function harmonicsson
figure; hold on;
for harmonic=1
    [artipuls,articheck]=artipulse(harmonic);
    arti=2350+articheck*150;
    arti=arti(512:length(arti));
    k=1;
    for i=1:length(arti)
        if k==4;
            g(round(i/4))=arti(i);
            k=1;
        end
        k=k+1;
    end
    pitch_data=jc_PitchData610son(artipuls,1024,1020,0.2,harmonic*2100,harmonic*2600,1); 
    plot(pitch_data/harmonic,'Color','black');
end
for harmonic=2
    [artipuls,articheck]=artipulse(harmonic);
    arti=2350+articheck*150;
    arti=arti(512:length(arti));
    k=1;
    for i=1:length(arti)
        if k==4;
            g(round(i/4))=arti(i);
            k=1;
        end
        k=k+1;
    end
    pitch_data=jc_PitchData610son(artipuls,1024,1020,0.2,harmonic*2100,harmonic*2600,1); 
    plot(pitch_data/harmonic,'Color','blue');
end
for harmonic=3
    [artipuls,articheck]=artipulse(harmonic);
    arti=2350+articheck*150;
    arti=arti(512:length(arti));
    k=1;
    for i=1:length(arti)
        if k==4;
            g(round(i/4))=arti(i);
            k=1;
        end
        k=k+1;
    end
    pitch_data=jc_PitchData610son(artipuls,1024,1020,0.2,harmonic*2100,harmonic*2600,1); 
    plot(pitch_data/harmonic,'Color','green');
end

plot(g,'Color','red');