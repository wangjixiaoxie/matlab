function edit_note

%edits notes or interval selected by mouse

%global variables for display
global h_main_amp onsets offsets labels

subplot(h_main_amp)
title('Add note: select left and right borders with mouse/ret or d to delete; q to quit','color',[1, 0, 0])

[left, right, button, flag]= get_xrange;

%if no change returned by get_xrange then exit
if flag == 0;
  title('')
  return
end

num_offs_pre = nnz(offsets < left);
num_ons_post = nnz(onsets > right);
num_ons_pre = nnz(onsets < left);
num_offs_post = nnz(offsets > right);

if  strcmp(char(button),'d')
   %if delete command eliminate notes/insert interval
   %copy leading unaffected notes
   for i = 1:num_ons_pre
     temp_ons(i,1)=onsets(i);
     temp_offs(i,1)= min(offsets(i),left);
     temp_labs(i,1)=labels(i);
   end
   
   %for the rest:
   for i = 1:num_offs_post
     index=length(onsets) - num_offs_post + i;
     t_index=num_ons_pre+i;  
     temp_ons(t_index,1)=max(right, onsets(index));
     temp_offs(t_index,1)=offsets(index);
     temp_labs(t_index,1)=labels(index);
   end
         
else
   %add/join/modify notes
   %copy leading unaffected notes
   for i = 1:num_offs_pre
     temp_ons(i,1)=onsets(i);
     temp_offs(i,1)=offsets(i);
     temp_labs(i,1)=labels(i);
   end
 
   %add edited note (which may encompass several notes)
   temp_ons(num_offs_pre + 1,1) = left;
   temp_offs(num_offs_pre + 1,1)= right;
   temp_labs(num_offs_pre + 1,1)= '0';

   %copy trailing unaffected notes
   for i = 1:num_ons_post
     index=length(onsets)-num_ons_post + i;
     t_index=num_offs_pre+1+i;  
     temp_ons(t_index,1)=onsets(index);
     temp_offs(t_index,1)=offsets(index);
     temp_labs(t_index,1)=labels(index);
   end

   %keep label if note # is unaltered
   if length(temp_ons) == length(onsets)
     temp_labs(num_offs_pre+1) = labels(num_offs_pre+1);
   end
end

onsets=temp_ons;
offsets=temp_offs;
labels=temp_labs';
  
subplot(h_main_amp)
title('')
disp_song;

