    %Calculates match rate
    check_stuff.(sofinterest).match.sum = sum(check_stuff.(sofinterest).match.values);
    check_stuff.(sofinterest).match.match_rate = [num2str(check_stuff.(sofinterest).match.sum(1)...
        ./check_stuff.(sofinterest).match.sum(2)*100) '%'];
    
    %Calculates timing of template matching
    check_stuff.(sofinterest).match.toff = [];
    for jj = 1:length(check_stuff.(sofinterest).match.trigger)
        check_stuff.(sofinterest).match.toff = [check_stuff.(sofinterest).match.toff; check_stuff.(sofinterest).match.trigger(jj).toffset];
    end
    
    check_stuff.(sofinterest).match.toff_mean = mean(check_stuff.(sofinterest).match.toff);
    check_stuff.(sofinterest).match.toff_sd = std(check_stuff.(sofinterest).match.toff);
    check_stuff.(sofinterest).match.toff_median = median(check_stuff.(sofinterest).match.toff);
    check_stuff.(sofinterest).match.toff_iqr = iqr(check_stuff.(sofinterest).match.toff);
    
    
    %Plots trigger timing vs. trigger number so you can see outliers
    lt_figure
    lt_plot(check_stuff.(sofinterest).match.toff);
    title(['Timing offsets for ' parameters.sofinterest ' (Offline match; correct label)'])
    xlabel('Syllable number')
    ylabel('Timing offset (ms)')
    saveas(figure(gcf), [SaveDir3 '/Toff_distro_' parameters.sofinterest '_(all)_' parameters.today_date], 'fig')
    
