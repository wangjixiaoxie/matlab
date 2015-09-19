function jc_plotabunch505(arrayfileD)
figure; hold on
spacing=180;
sample_min=1;
sample_max=25;
ymin=2000;
bulge=400; % Fits stuff on graph
ymax=8000; %ymin+bulge+spacing*(sample_max-sample_min);

for i=sample_min:sample_max
    plotshiftedD=arrayfileD(i).pitches(1:200)+spacing*i;
    plot(plotshiftedD); xlim([0 200]); ylim([ymin ymax]); title('Directed')
end
