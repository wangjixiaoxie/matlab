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
		if ~feof(fin)
		[fn, count] = fscanf(fin,'%s',1);
		else
		break
		end
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

sums = sum(tran_mat,2);
for ii = 1:length(sums)
	if (sums(ii) ~= 0)
		tran_mat(ii,:) = 100.0*tran_mat(ii,:)./sums(ii);
	end
end

reps = cell([length(notes),1]);
fin = fopen(batch,'r');
while (~feof(fin))
	fn = fscanf(fin,'%s',1);
	prev_note='+';prev_ind=1;
	if (exist(fn,'file'))
		load(fn);
%		disp(fn);
		for ii = 1:length(notes)
			nt=notes(ii);
			labels = lower(labels);
			p = find(labels==nt);
			if (length(p)>0)
				tmp_reps = [];
				cnt=1;
				for kk = 2:length(p)
					if (p(kk)==p(kk-1)+1)
						cnt = cnt + 1;
				        else
						tmp_reps = [tmp_reps,cnt];
						cnt=1;
				        end
			        end
				tmp_reps=[tmp_reps,cnt];

				%load up the cell array
				reps1 = reps{ii};
				if (length(reps1)<max(tmp_reps))
					vtmp = reps1;
					reps1 = zeros([1,max(tmp_reps)]);
					reps1(1:length(vtmp)) = vtmp;
			        end
				for jj = 1:length(tmp_reps)
				   reps1(tmp_reps(jj))=reps1(tmp_reps(jj))+1;
			        end
				reps{ii} = reps1;
			end
		end
	end
end
fclose(fin);
return;
