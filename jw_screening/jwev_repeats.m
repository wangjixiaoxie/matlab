function [tran_mat,reps,notes,nums]=jwev_repeats(batch,ALL_LOWER);
%[tran_mat,reps,notes,nums]=jwev_repeats(batch,ALL_LOWER);
% analyzes data
% if ALL_LOWER=1 then will consider all notes in .not.mat files to be lower case
%diff from ev_repeats in that rep trains are treated as singlets.

% jw commented out display of filenames (2 changes of disp(fn))
if (~exist('ALL_LOWER'))
	ALL_LOWER = 0;
end

notes = ['+/'];
tran_mat = zeros(length(notes));
nums=[0,0];
%for ignoring one file listed in the batch randomly each time
ignore=rand;
batchlength=load_batchf(batch);
ignore=round(ignore*length(batchlength));
fcounter=0;

fin = fopen(batch,'r');
while (~feof(fin))
	fcounter=fcounter+1;
	fn = fscanf(fin,'%s',1);
	if fcounter==ignore
		continue;
		%the old code.. delete?
		%if ~feof(fin)
		%[fn, count] = fscanf(fin,'%s',1);
		%else
		%break
		%end
	end
	if (exist(fn,'file'))
%		disp(fn);
		load(fn);
	%change rep trains to singlets
	todel=[];
	for rs=2:length(labels)
		if labels(rs)==labels(rs-1)
			todel=[todel rs];
		end	
	end
	labels(todel)=[];
        %exclude notes labelled as 0
        labels(findstr(labels,'0'))=[];
		prev_note = '+';
		prev_ind  = 1;
		nums(1) = nums(1) + 1;
		for ii = 1:length(labels)
			if (ALL_LOWER==1)
				nt = lower(labels(ii));
			else
               nt = labels(ii);
            end
				
			pos = find(notes==nt);
			if (length(pos)==0)
				notes = [notes,nt];
				temp = tran_mat;
				ln = size(tran_mat,1);
				tran_mat = zeros(ln+1);
				tran_mat(1:ln,1:ln) = temp;
				nums = [nums,0];
				pos = ln+1;
			end
			nums(pos) = nums(pos) + 1;
			tran_mat(prev_ind,pos) = tran_mat(prev_ind,pos) + 1; 
			prev_ind = pos;
		end
		pos = find(notes=='/');
		nums(pos) = nums(pos) + 1;
		tran_mat(prev_ind,pos) = tran_mat(prev_ind,pos) + 1; 
	end
end
fclose(fin);

sums = sum(tran_mat,2) %sums should equal nums??
nums
for ii = 1:length(sums)
	if (sums(ii) ~= 0)
		tran_mat(ii,:) = 100.0*tran_mat(ii,:)./sums(ii);
	end
end
