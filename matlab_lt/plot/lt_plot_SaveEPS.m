SaveDir='/bluejay4/lucas/Dropbox/SCIENCE/BRAINARD_LAB/MANUSCRIPTS/SeqDepPitch/FIGURES/MATLAB';
timestamp=lt_get_timestamp(0);
randnum=randi(10000, 1);
SavePath=[SaveDir '/' timestamp '_rand' num2str(randnum) '.eps'];


print('-depsc2', '-noui', '-adobecset', '-painters', SavePath)

%% Using export_fig, poptentialyl better?

% -transparent

if (0) % not working so far. might need to download ghostscript
savetemp=[SaveDir '/export_fig_version.eps'];


export_fig('-painters', '-depsc', '-eps', '-CMYK', savetemp);
export_fig savetemp -png
end