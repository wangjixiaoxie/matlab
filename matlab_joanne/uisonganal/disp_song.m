 function [h_amp_plot, h_labels] = disp_song 
    %displays amplitude and label data in current plot window
    %displays data in the range set by axis limits saved with current axis handle
    % Now needs handles of amplitude and label axes since the are
    % separate now (BDW)
    
    %global variables for display
    
    global onsets offsets labels 
    global d_song time_d_song
    global h_main_amp h_main_labels
    
    %get value for scaling note display to song amplitude
    ymax=max(d_song);
    
    %display     
    set(gcf, 'CurrentAxes', h_main_amp)
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

    %plot labels
    
    set(gcf, 'CurrentAxes', h_main_labels)
    % Clear axes and then display new text 
    h_old_labels = findobj(h_main_labels,'Tag','SegLabels');
    if ( ~isempty(h_old_labels))
      delete(h_old_labels);
    end
    h_labels = text(onsets,ones(size(onsets))*0.5,labels','Clipping','on');
    set(h_labels,'Tag','SegLabels')
    
    set(gcf, 'CurrentAxes', h_main_amp)

    
