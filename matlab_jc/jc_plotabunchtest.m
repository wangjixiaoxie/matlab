function [medpitch,boomed]=jc_plotabunchtest(arrayfileD,f,xmax)

%notelength=850;
%Parameters
spacing=500;
note_length=6400;

sampling_rate=32000;

ymin=0;
ymax=20000; 

xmin=0;
%xmax=length(arrayfileD(1).pitches);

%Smooth the data
for i=1:length(arrayfileD)
    [holder]=SmoothData(f(i).datt,sampling_rate,1);
    smooth(i).smoothed=holder;
end

%Get average smoothed note

figure; hold on
for i=1:length(arrayfileD)
    processed=jc_pinterp(arrayfileD(i).pitches);
    medpitch(i)=median(processed(500:600));
    boomed(i)=std(processed(350:450));
    processed=processed+spacing*i;
    plot(processed); xlim([xmin xmax]); ylim([ymin ymax]); title('1-10')
    plotted(i).trace=processed-380*i;
end
