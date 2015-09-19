%% LT 4/6/15 - 
% e.g.:
% batch = 'batch.keep';
% sort_method = 'date' (sorts by time/date);

function lt_sort_batch(batch,sort_method)

disp('NOTE: this currently works only for evtafv4 filenames, using date, and only actually looks at time')

fid=fopen(batch);
tmp=textscan(fid,'%s');

song_list=tmp{1};

% convert names to times
datenumlist=[];
for i=1:length(song_list);
    
    % find index holding date and time
    undscore=strfind(song_list{i},'_');
    
    datetime=song_list{i}(undscore(end)-6:undscore(end)+6);

    % put all dates into a list
    datenumlist(i)=datenum(datetime,'ddmmyy_HHMMSS');
    
end


% sort that list of dates
[~, Inds]=sort(datenumlist);

% reorder song filenames based on date
song_list_ordered=song_list(Inds);

% put back into new file
batch2=[batch '.sorted_' sort_method];
fid2=fopen(batch2,'w');

for i=1:length(song_list_ordered);
    fprintf(fid2,'%s\n',song_list_ordered{i});
end

fclose(fid2);

end



