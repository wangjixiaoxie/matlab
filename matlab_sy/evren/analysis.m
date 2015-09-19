function [tran_mat,reps,notes,nums]=analysis(batch,ADDNOTMAT);
% [tran_mat,reps,notes,nums]=analysis(batch,ADDNOTMAT);
% analyzes data

if (~exist('ADDNOTMAT'))
    ADDNOTMAT=0;
end

notes = ['+/'];
tran_mat = zeros(length(notes));

nums=[0,0];
fin = fopen(batch,'r');
while (1)
	fn = fgetl(fin);
	if (~ischar(fn))
		break;
    end
    if (ADDNOTMAT==1)
        fn = [fn,'.not.mat'];
    end
	if (exist(fn,'file'))
		%disp(fn);
		load(fn);
		prev_note = '+';
		prev_ind  = 1;
		nums(1) = nums(1) + 1;
		for ii = 1:length(labels)
			nt = lower(labels(ii));
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
while (1)
	fn = fgetl(fin);
	if (~ischar(fn))
		break;
	end
	prev_note='+';prev_ind=1;
    if (ADDNOTMAT==1)
        fn = [fn,'.not.mat'];
    end
	if (exist(fn,'file'))
		load(fn);
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
