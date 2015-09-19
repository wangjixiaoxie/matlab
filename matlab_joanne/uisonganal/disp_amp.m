 function h_amp_plot = disp_amp 
    %displays amplitude and label data in current plot window
    %displays data in the range set by axis limits saved with current axis handle
    %global variables for display
    
    global onsets offsets labels 
    global d_song time_d_song
    
    %get value for scaling note display to song amplitude
    ymax=max(d_song);
    
    %display     
    cla
    hold on
    
    %plot song
    h_amp_plot = plot(time_d_song, d_song);
    set(h_amp_plot,'Tag','AmpPlot');
    
    %plot notes
    plot([onsets'; onsets'],[ymax*ones(size(onsets'));zeros(size(onsets'))],'m');
    plot([offsets'; offsets'],[ymax*ones(size(offsets'));zeros(size(offsets'))],'m');
    plot([onsets';offsets'],[ymax*ones(size(onsets'));ymax*ones(size(onsets'))],'m');

    hold off

    
