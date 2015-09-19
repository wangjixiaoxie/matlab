function jc_plotabunch501(arrayfileD,arrayfileUD,arrayfileL)
figure; hold on
subplot(131); hold on
spacing=100;
sample_min=1;
sample_max=17;
ymin=2600;
bulge=400; % Fits stuff on graph
ymax=ymin+bulge+spacing*(sample_max-sample_min);

for i=sample_min:sample_max
    plotshiftedD=arrayfileD(i).pitches(1:170)+100*i;
    plot(plotshiftedD); xlim([0 200]); ylim([ymin ymax]); title('Directed')
end
subplot(132); hold on
for i=sample_min:sample_max
    plotshiftedUD=arrayfileUD(i).pitches(1:170)+100*i;
    plot(plotshiftedUD); xlim([0 200]); ylim([ymin ymax]); title('Undirected')
end
subplot(133); hold on
for i=sample_min:sample_max
    plotshiftedL=arrayfileL(i).pitches(1:170)+100*i;
    plot(plotshiftedL); xlim([0 200]); ylim([ymin ymax]); title('Lesioned')
end