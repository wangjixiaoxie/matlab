function save_data(fname, pathname, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win)


%save format with data
%format_index=findstr('.not1',fname);
%if isempty(format_index)
%   format=0;
%else
%   format=1;
%end

note_file_full  = fullfile(pathname,fname);
disp(['writing note data to ', note_file_full])
%eval(['cd ', pathname]);
save(note_file_full, 'Fs','onsets','offsets','labels',...
    'threshold','min_int','min_dur','sm_win');


