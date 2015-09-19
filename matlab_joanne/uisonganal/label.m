function label(onsets,offsets)

%allows modification of note labels

%global variables for display
global d_song time_d_song h_main_amp h_main_spect h_main_labels labels
global win_percent

disp('Use left mouse button to select note for labelling')
disp('then type labels with pointer in plot window')
disp('labels must be single character')
disp('return skips to next note to right')
disp('backspace skips to next note to left')
disp('type q to quit')

%by default, start with first note
note_num=1;
subplot(h_main_amp)
[h_amp_plot, h_labels]=disp_song;
title('LABEL MODE! (q to quit)','color',[1,1,0]);      
set(h_labels(note_num),'color',[1 0 0]);

%change_flag keeps track of whether notes have been typed since last return
change_flag = 0;

subplot(h_main_amp);

while 1
   
   %get input from keyboard or mouse
   [t,y,key]=ginput(1);

   if ~isempty(key) & key == 1
     xlim = get(gca,'xlim'); 
     if t < xlim(1)
       move_left(win_percent)
       subplot(h_main_spect);
       move_left(win_percent)
       subplot(h_main_labels);
       move_left(win_percent)
       subplot(h_main_amp);            
     elseif t > xlim(2)
       move_right(win_percent)
       subplot(h_main_spect);
       move_right(win_percent)
       subplot(h_main_labels);
       move_right(win_percent)
       subplot(h_main_amp);
     else
       note_num=max(nnz(t>onsets),1);
       %display  labels & current insertion point
       subplot(h_main_amp)
       [h_amp_plot,h_labels]=disp_song;
       title('LABEL MODE! (q to quit)','color',[1,1,0])      
       set(h_labels(note_num),'color',[1 0 0]);
     end 
   elseif strcmp(char(key),'q')
     subplot(h_main_amp)
     [h_amp_plot,h_labels]=disp_song;
     title('')          
     break
   elseif isempty(key) | key == 13  			%return
     if change_flag == 1
       change_flag = 0;
     else
       %advance insertion point
       note_num=min(note_num+1,length(labels));
     end  
     %display  new labels & current insertion point
     subplot(h_main_amp)
     [h_amp_plot,h_labels]=disp_song;
     title('LABEL MODE! (q to quit)','color',[1, 1, 0])      
     set(h_labels(note_num),'color',[1 0 0]);        
   elseif key == 8 | key == 127    %backspace or delete
     note_num=max(1,note_num-1);
     subplot(h_main_amp)
     [h_amp_plot,h_labels]=disp_song;
     title('LABEL MODE! (q to quit)','color',[1, 1, 0])      
     set(h_labels(note_num),'color',[1 0 0]);         
   else       
     labels(note_num)=char(key);
     note_num=min(note_num+1,length(labels));
     change_flag=1;       
   end
 end
