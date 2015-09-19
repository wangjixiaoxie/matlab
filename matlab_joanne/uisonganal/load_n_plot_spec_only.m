function load_n_plot_spec_only(filt)

% function load_n_plot_spec_only(filt)
% 
% This function calles loadraw and plotspec. FILT is a filter
% string for use by loadraw.  

[song,fs,name] = loadraw(filt);

plotspec(song,fs,name)