% label songs
% change miss label
% look at template of miss and hit
% find most extreme syllables in missed ones, and use that to make a new template. 
% you want to be able to hit the ones that you would have missed before.  so change AND NOT to OR YES logic for the old targets.
% 
% look at individual slices for:
% 1) hits, misses.  divide the latter into those you wanted to miss, andt hos not.  make template from fringe of what you wanted, non-fringe, and those you wanted. 
% if labeled both catch and notcatch, combine them by first changing misslables to same
%     thing, and then using transfer calls.
% then use lt_count, since all labeled songs are now in batch.catch.keep.
% 
% make new template
% load old one for comparison.
% 
% save template and modify old template and parameters.
% 
% test with uievtafsim and lt_db_check...


%% making new template
%load old one
load('../pu13bk43b_harmonics_3.dat');

%make new one
pu13bk43b_harmonics_4=pu13bk43b_harmonics_3 ;
pu13bk43b_harmonics_4(:,[8:10])=[];

pu13bk43b_harmonics_4=[pu13bk43b_harmonics_4 template_temp pu13bk43b_harmonics_3(:,[8:10])];

%save template
wrt_templ('pu13bk43b_harmonics_4.dat', pu13bk43b_harmonics_4)

%move to root folder
movefile('pu13bk43b_harmonics_4.dat','..')

%modify the template_parameters structure (which is used in
%lt_db_check_template_timing_and_hit_rate)
lt_add_to_template_parameters







