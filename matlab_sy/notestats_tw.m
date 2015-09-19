function [b_labels] = notestats_tw(bt)

%[b_labels, b_onsets, b_offsets] = notestats
%
%   reads a series of notefiles
%   returns the number of distinct notes,
%        a matrix containng note durations, one row per disitnct note
%        a string representing all song labels


%until the end of file, load the fiels and concatenate onsets offsets and labels
% into  single column vectors with '/' separating each input vector in label column


%get name of metafile containing notefile names


   meta_fid=fopen(bt);
   if meta_fid == -1 
      disp('cannot open file' );
      disp (metafile);
   end


%the special symbol '/' is used to separate songs

b_onsets=[0];
b_offsets=[0];
b_labels=['/'];


while 1
   %get notefile name
     notefile = fscanf(meta_fid,'%s',1);
     
     outnotefile=[notefile '.not.mat']
   %end when there are no more notefiles 
     if isempty(notefile);
        disp('End of notefiles');
        break
     end
   
%    if notefile exists, get it
     if exist(outnotefile);   
       load(outnotefile);
     else
       disp(['cannot find ', outnotefile]);
       labels=[];
     end
    
   %make sure labels is a column vector
   if(~isempty(labels))
        [x,y]=size(labels);
        if y > x
            labels=labels';
        end 
   
   %append input strings to batch strings and terminate with '/'
     b_onsets = [b_onsets; onsets; 0];
     b_offsets = [b_offsets; offsets; 0];
     b_labels = [b_labels; labels; '/'];
   end
end
     
% %find all the unique notes and their frequency of occurence
% [elts, nelts] = elements(b_labels);
% 
% %return elts as a column vector
% elts = elts';
% 
% 
% 
% 
% %return 
% 
% for i = 1:length(elts);
%    cur_durs = get_durs(b_onsets, b_offsets, b_labels, elts(i));
%    cur_durs = cur_durs';
%    dur_data(i,1:length(cur_durs)) = cur_durs;
% end
% 
% 
% %plot a summary of note durations
% plot_notes(elts,elts,dur_data);xlabel('note');ylabel('duration');
% end