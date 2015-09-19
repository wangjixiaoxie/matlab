function [lin,scon,st,lcon,entr]=jw_stereot5(batch,condition,rep)
% [lin,scon,st,lcon,entr]=jw_stereot5(batch,condition,rep)
%
% The fully-functioning version of the jw_stereot code.
%
% Can treat repeat trains as single syls
% or exactly as sung, with the rep argument.
%
%For bootstrapping to check for influence of single files,
%one input file from the batchfile can be omitted randomly each time.
%
%Included some ev_repeats code so should be stand-alone.
%
% batch should have one filename per line only. Each file should
% contain a matlab variable called "labels," which is a vector of
% characters representing syllables in the order sung. Any syllables
% labeled '0' will be excluded. Multiple song beginnings and endings
% can be notated within one 'labels' variable with '+' and '/'
% characters respectively.
%
% condition = integer 1 or 2.	(1)=include endings
%				(2)=disregard endings ('internal' values only)
%                                   (see references below)
%   'disregard endings' means the variability in which syls can end song does
%   not affect stereotypy computations ('to-end' transitions discarded).
%
% rep = 's' or 'r' to indicate if repeat trains should be collapsed
% so as to be treated as a single syllable ('s') or unchanged as a
% train of repeats ('r')
%
% Many stereotypy values are computed and reported in standard output,
% first per-syllable and then whole-song, but only whole-song values
% that fit the condition argument are returned by the function.
%
% First, per-syllable values are reported, wherein only the values
% pertaining to the 'condition' argument are displayed.
% LC is reported PER-NOTE as LC value, not normalized. 
% Entropy is reported PER-NOTE as entropy (p*log2(p)), not
% normalized.
%
% Finally, whole-song values are shown. LC and entropy are
% prevalence-weighted sums. Values are reported with or without
% 'to-end' transitions, so you can compare. The subset of these
% displayed values which matches the 'condition' argument
% is returned by the function.
%
% lin = either 'sequence linearity' (with endings) or internal linearity
%	(excludes endings) according to 'condition' arg.
%	see Scharff & Nottebohm 1991, Foster & Bottjer 2001
% scon= 'sequence consistancy.'  see Scharff & Nottebohm 1991.
% st  = average of lin and scon, see Scharff & Nottebohm 1991.
% lcon= JW's 'linear consistancy' measure. scores each syllable as
%	having a linear consistancy of 0 - 1. If a syl has 2 transitions
%	with .99 and .01 probabilities, lc is almost 1,
%	and increasing either the prevalence of the minor trans or
%	the # possible trans will smoothly decrease the score.
%	lc = 1 - 1/(sum(p/(1-p))). Whole-song lcon is weighted by syl prevalence.
%       No reference, shown for comparison only.
% entr= entropy. Weighted sum of entropy (p*log2(p)) for each syllable.
%	Whole-song value is weighted by syllable prevalence.
%

rand('state',sum(100*clock)); %for bootstrapping-- currently set up to omit one file each time. Note the "ignore=0" line below disables this feature. Must comment that out to use.
if (~exist('INTRO'))
	INTRO='+';
end

%Changing condition number to match modified code.. 
if condition>2; return; end;
if condition==2
	condition=3;
end

%%%%%%%%%%%%%%%%%%%%%%%
% These lines from ev_repeats.m are included to replace a call to that function.
% [tmat,reps,nt,nm]=jwev_repeats(batch,1);

notes = ['+/'];
tran_mat = zeros(length(notes));
nums=[0,0];
%for ignoring one file listed in the batch randomly each time
ignore=rand;
batchlength=load_batchf(batch);
ignore=round(ignore*length(batchlength));
fcounter=0;

%this disables the 1-file exclusion for testing. Comment out to exclude.
ignore=0;
%

fin = fopen(batch,'r');
while (~feof(fin))
	fn = fscanf(fin,'%s',1);
	fcounter=fcounter+1;
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
	if rep=='s'
	%change rep trains to singlets.
	todel=[];
		for rs=2:length(labels)
			if labels(rs)==labels(rs-1)
				todel=[todel rs];
			end	
		end
		labels(todel)=[];
	end
        %exclude notes labelled as 0
        labels(findstr(labels,'0'))=[];
		prev_note = '+';
		prev_ind  = 1;
		nums(1) = nums(1) + 1;
		for ii = 1:length(labels)
		% All characters in 'labels' variable will be converted to their lowercase
		% counterparts. Un-comment here to build back in the ability to treat upper and
		% lower separately.
			%if (ALL_LOWER==1)
				nt = lower(labels(ii));
			%else
               		%	nt = labels(ii);
            		%end
			
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
			if nt=='+' %If labeled song includes a new start, automatically include an end
				   %and don't log anything as transitioning to this '+'
			nums(find(notes=='/'))=nums(find(notes=='/'))+1;
			prev_ind = pos;
			continue;
			end
			tran_mat(prev_ind,pos) = tran_mat(prev_ind,pos) + 1;
			if nt=='/'
			prev_ind = 1; %if a '/' (song end marker) was included in a song, automatically include a start
			nums(1) = nums(1) + 1;
			continue
			end 
			prev_ind = pos;
		end
		if labels(ii)~='/'
		pos = find(notes=='/');  
		nums(pos) = nums(pos) + 1;
		tran_mat(prev_ind,pos) = tran_mat(prev_ind,pos) + 1;
		end
	end
end
fclose(fin);

tsums = sum(tran_mat,2); %tsums should equal nums, except for '/' label which by def is the end and
                       %therefore transitions to nothing.
			%so tsums not nums should be used for converting tmat back and forth b/w
			%probabilities and raw counts.
for ii = 1:length(tsums)
	if (tsums(ii) ~= 0)
		tran_mat(ii,:) = 100.0*tran_mat(ii,:)./tsums(ii);
	end
end

%%%%%%%%%%%%%%%%%%%%%%% end code from jwev_repeats.m
tmat=tran_mat;
nt=notes;
nm=nums;

outs=[0 0 0 0 ; 0 0 0 0;0 0 0 0;0 0 0 0];	%rows for linearity, consistancy, lc, and entr. first column with endings, 3nd = without. 2nd and 4th are archaic--ignore.

% b= # files listed in batch file, minus the one ignored file
if ignore ~=0
b=length(batchlength)-1;
else
b=length(batchlength);
end
fprintf(1,['\n----------Stereotypy: ',batch,'------------\n\r',num2str(b),' songs to analyze.\n\r']);


% Express tmat as percentages, and get raw numbers tmat. Stereot4 may need this code.
ptmat=tmat.*.01;
for z=1:size(ptmat,1)
raw_tmat(z,:)=ptmat(z,:).*tsums(z); %use tsums instead?
end
sylnums=[];
for ii=1:size(tmat,1);
	if (~strcmp(nt(ii),'0'))&(~strcmp(nt(ii),'+'))&(~strcmp(nt(ii),'/'));
		sylnums=[sylnums,nm(ii)]; %new ev_repeats-style 'nums' vector for only those syllables meeting conditions of this IF loop
	end
end

%%%%
% calculations begin here...
total_domtrans=[0,0,0,0];lin_con_sum=[0 0 0 0];entr_sums=[0 0 0 0];  % The 4 values are from previous code.. here we use the 1st and 3rd only. 1st=w/endings 3rd=without.
for ii=1:size(ptmat,1);
	        if (~strcmp(nt(ii),'0'))&(~strcmp(nt(ii),'+'))&(~strcmp(nt(ii),'/'));
			                %(~strcmp(nt(ii),'i'))&(~strcmp(nt(ii),'j')
					                total_domtrans=[(total_domtrans(1)+max(raw_tmat(ii,:))),0,(total_domtrans(3)+max(raw_tmat(ii,3:size(raw_tmat,2)))),0];       %sums of # dom transitions sung: [withendings,ignore,without endings,ignore]

% First, LC and entr. without endings and with.
		for g = 1:2
			if g==2;g=3;end;  %this code only uses 2 options for "condition," the original 1 and 3 conditions.
			sum_thissyl_lincon=0.0;
			sum_thissyl_entr=0.0;
			if g == 1
			mat=ptmat; s=1; cond='with-endings'; rmat=raw_tmat;
			end
			if g== 3
			mat=ptmat; s=3;cond='exclude-ends (internal transitions only)'; rmat=raw_tmat;
			end
			todo=find(mat(ii,s:length(mat(ii,:))));
			if length(todo)<=1
				sum_thissyl_lincon=1;
			% 	so perfectly stereotyped transitions get a '1' for LC and keep the '0' for entropy
			else
				for inrow=1:length(todo);
					sum_thissyl_lincon=sum_thissyl_lincon+(mat(ii,todo(inrow)+(s-1))/(sum(mat(ii,s:length(mat(ii,:))))-mat(ii,todo(inrow)+(s-1))));
					sum_thissyl_entr=sum_thissyl_entr-(mat(ii,todo(inrow)+(s-1))*(log2(mat(ii,todo(inrow)+(s-1)))));
					% the "s-1" term is because todo does not hold absolute references to column positions in ptmat, it is references to the result of a find command on a vector with a starting index of 1 or 3. (so that condition 1 or 2 could have different find commands) Should have just always started at 1 and then edit the resulting vector according to the 'condition' being evaluated  
				end
			sum_thissyl_lincon=1-1/sum_thissyl_lincon;
				
			end
			%prevalence weight the values found for this syl, and you have the contribution to the whole-song score.
			outs(3,g)=outs(3,g) + (sum(rmat(ii,s:size(rmat,2)))/sum(sum(rmat(3:size(rmat,1),s:size(rmat,2)))))*(sum_thissyl_lincon);
			outs(4,g)=outs(4,g)+(sum_thissyl_entr*(sum(rmat(ii,s:length(rmat(ii,:))))/sum(sum(rmat(3:size(rmat,1),s:size(rmat,2))))));
			%outs(3,g)=outs(3,g) + (nm(ii)/sum(sylnums))*(sum_thissyl_lincon);
			%outs(4,g)=outs(4,g)+(sum_thissyl_entr*(nm(ii)/sum(sylnums)));
			if g==condition
				disp_cond=cond; % for display at the end of output
				fprintf(1,['Syllable ''',nt(ii),''': condition = ',cond,'.\n\t\b   total transition types = ',num2str(length(find(mat(ii,s:length(mat(ii,:))))))])
				if s==1
			%	if g==1
				if mat(ii,2)~=0
					fprintf(1,[' (*''to-end''*; ',num2str(100*((rmat(ii,2))/sum(rmat(:,2)))),'%% of ends.)'])
					%flags syllables that are song-enders, tells you percent of all endings that this syllable accounts for.
				end
				end
				clear maxval;clear c1;
				[maxval,c1]=max((mat(ii,:)));
				fprintf(1,['\n             Dominant transition = ',num2str(100*maxval),'%%, to ''',num2str(nt(c1)),'.''\n            Transitional entropy = ',num2str(sum_thissyl_entr),'\n              Linear Consistancy = ',num2str(sum_thissyl_lincon),'\n            Prevalence = ',num2str(sum(rmat(ii,s:size(rmat,2)))/sum(sum(rmat(3:size(rmat,1),s:size(rmat,2))))),'\n\n']);
			end
				%sum_thissyl_entr %comment out! for debug only
				%sum_thissyl_lincon %comment out! for debug only
				
		end
	end
end
fprintf(['  Syllables included for analysis: ',num2str(length(sylnums)),'.\n  Should match total song syls labeled: ',num2str(length(nt)-2),' syllables.\n']);
fprintf(['\nDominant transition types\n       trans types sung:\t',num2str(total_domtrans(1)),'\nother than to-end types:\t',num2str(total_domtrans(3)),',\n\nTotal trans types\n       trans types sung:\t',num2str(sum(sum(raw_tmat(3:size(raw_tmat,1),:)))),'\nother than to-end types:\t',num2str(sum(sum(raw_tmat(3:size(raw_tmat,1),3:size(raw_tmat,2))))),'\n\n'])

%total syl types minus 1, divided by total trans NOT FROM begin or end... = il
outs(1,3)=(size(ptmat,1)-3)/length(find(ptmat(3:size(ptmat,1),3:size(ptmat,2))));
% traditional sl for reference
outs(1,1)=(size(ptmat,1)-2)/length(find(tmat(3:size(tmat,1),:)));
% traditional sc
outs(2,1)=total_domtrans(1)/sum(sum(raw_tmat(3:size(raw_tmat,1),:)));
% sc disregarding to-end transitions-- "iSC"
outs(2,3)=total_domtrans(3)/sum(sum(raw_tmat(3:size(raw_tmat,1),3:size(raw_tmat,2))));
% outs(3,:) already defined in 'for' loop. JW-style Linear Consistancy, weighted by syllable prevalence. To see unweighted, go to last for loop and swap commented line for real line when defining lin_con_sum
% outs(4,:) already defined in 'for' loop. Unnormalized entropy, weighted by syllable prevalence
%ptmat,loopy_ptmat
%raw_tmat, loopy_tmat
fprintf(1,['\n\t  Whole-song values\n\n     SL = ',num2str(outs(1,1)),'\n     IL = ',num2str(outs(1,3)),'\n\n     SC = ',num2str(outs(2,1)),'\n    iSC = ',num2str(outs(2,3)),'\n\n     LC = ',num2str(outs(3,1)),'\n    iLC = ',num2str(outs(3,3)),'\n\n   entr = ',num2str(outs(4,1)),'\n i-entr = ',num2str(outs(4,3)),'\n'])
%disp(['For the ''',disp_cond,''' condition:'])
lin=outs(1,condition);scon=outs(2,condition);st=(lin+scon)/2;lcon=outs(3,condition);entr=outs(4,condition);
fprintf(1,['\nYou chose the ''',disp_cond,''' condition:\n','\n          Lin = ',num2str(lin),'\n       Consis = ',num2str(scon),'\n their avg=ST = ',num2str(st),'\n        JW_LC = ',num2str(lcon),'\n  Entropy sum = ',num2str(entr),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
%fprintf(1,['\nFor the ''',disp_cond,''' condition:\n','\n       JW_LC = ',num2str(lcon),'\n          ST = ',num2str(st),'\n Entropy sum = ',num2str(entr),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
%disp(['          SL = ',num2str((outs(1,condition)+outs(2,condition))/2)])
%fprintf(1,[Entropy sum = ',num2str(entr_sums(condition)),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
