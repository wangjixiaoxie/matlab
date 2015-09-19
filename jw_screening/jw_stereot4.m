function [lin,scon,st,lcon,entr]=jw_stereot4(batch,condition,rep,INTRO)
% [lin,scon,st,lcon,entr]=jw_stereot4(batch,condition,rep,intro)
%
%unlike stereot2 and 3, this can treat repeat trains as single syls
%or exactly as sung, with the rep argument.
%For bootstrapping to check for influence of single-outliers,
%one input file from the batchfile is omitted randomly each time.
%Included some ev_repeats code so no other m-files needed.
%
% batch should have one filename per line only. Each file should
% contain a matlab variable called "labels," which is a vector of
% characters representing syllables in the order sung. Any syllables
% labeled '0' will be excluded. Multiple song beginnings and endings
% can be notated within one 'labels' variable with '+' and '/'
% characters respectively.
%
% condition = integer 1 thru 4.	(1)=withendings, nonloopback.
%				(2)=withendings,loopback-corrected.
%				3 and 4 are without endings
%   Combinations of 2 different features:
%   1.'without endings' means the variability in which syl ends a song does
%   not affect stereotypy computations ('to-end' transitions discarded).
%   2.loopback-correction tries to eliminate a decrement in stereotypy
%   for songs that have variable numbers of highly-stereotyped
%   motifs. eg if there are 3 files with the labels variables
%   abcd, abcd, and abcdabcd, loopback correction makes the reported
%   stereotypy values the same as for 3 files abcd, abcd, and abcd.
%   Thus repeating the motif does not decrement stereotypy, but
%   motif variants like abcabcd or abcdbcd do introduce a decrement.
%   Furthermore, those two variants introduce only a reduced decrement,
%   unlike variants like abcbcd, which are not affected at all by loopback
%   correction because they involve restarts that occur 1) after syllables
%   that never end song AND 2) before syllables that never start song.
%   (See descrip of "intro" argument, below.)
%
% rep = 's' or 'r' to indicate if repeat trains should be collapsed
% so as to be treated as a single syllable ('s') or unchanged as a
% train of repeats ('r')
%
% intro is a syllable label that represents an intronote, for the
% purpose of determining loopback transitions.
%   (default is the + label, which is added to each file's labels variable.
%   Loopback score for a transition depends both on the "ending propensity"
%   for the "from" syllable and on the "starting propensity" of the "to"
%   syllable of the transition in question. So "starting propensity"
%   can be judged based on either how often the song begins with a
%   certain label or how often the POST-INTRONOTE song begins with
%   that label. For loopback correction, songs are considered to
%   start with whatever label follows the "intro" argument.)
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
% 'to-end' transitions and with or without loopback-correction. The
% subset of these displayed values which matches the 'condition' argument
% is returned by the function.
%
% lin = 'sequence linearity' (with endings) or internal linearity
%	(excludes endings) according to 'condition' arg.
%	see Scharff & Nottebohm 1991, Foster & Bottjer 2000
% scon= 'sequence consistancy.'  see Scharff & Nottebohm 1991.
% st  = average of lin and scon, see Scharff & Nottebohm 1991.
% lcon= JW's 'linear consistancy' measure. scores each syllable as
%	having a linear consistancy of 0 - 1. If a syl has 2 transitions
%	with .99 and .01 probabilities, lc is almost 1,
%	and increasing either the prevalence of the minor trans or
%	the # possible trans will smoothly decrease the score.
%	lc = 1 - 1/(sum(p/(1-p))). Whole-song lcon is weighted by syl prevalence.
% entr= entropy. Weighted sum of entropy (p*log2(p)) for each syllable.
%	Whole-song value is weighted by syllable prevalence.
%

rand('state',sum(100*clock)); %for bootstrapping-- currently set up to omit one file each time.
if (~exist('INTRO'))
	INTRO='+';
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

%this disables the 1-file exclusion for testing
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
		labels					%%%%%delete%%%%
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

outs=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];	%rows for linearity, consistancy, lc, and entr.
					%columns = the 4 conditions

% b= # files listed in batch file, minus the one ignored file
b=length(batchlength)-1;
fprintf(1,['\n----------Stereotypy: ',batch,'------------\n\r',num2str(b),' songs to analyze.\n\r']);
if isempty(findstr(INTRO,nt))
	fprintf(1,['The requested intronote ''',nt(findstr(INTRO,nt)),''' does not occur.\nUsing the true starts of songs as loop-starters.\n'])
	INTRO='+';
end

%starters= row of tmat for 'intro' syllable, expressed as percents, with 100% being all NON-i-TO-i transitions only
% i.e. what does INTRO transition to? Will show which syls start the song (where start means what comes after INTRO)

starters=[];
starters=.01.*tmat(findstr(INTRO,nt),:).*nm(findstr(INTRO,nt));	%'raw starters'
starters(findstr(INTRO,nt))=0;
starters=starters./sum(starters);

% Begin to make loopback-correction scores matrix tmat, count good syllables..
ends=[];loopback1=[];sylnums=[];
for ii=1:size(tmat,1);
	% a non-percentage-based tmat:
	raw_tmat(ii,:)=.01.*tmat(ii,:).*nm(ii); %use sums instead?
	if (~strcmp(nt(ii),'0'))&(~strcmp(nt(ii),'+'))&(~strcmp(nt(ii),'/'));
			%(~strcmp(nt(ii),'i'))&(~strcmp(nt(ii),'j')
		ends=[ends,(raw_tmat(ii,2))]; %how many song endings does each syllable do? sum is total # of endings..
		thissyl_trans_without=find(tmat(ii,3:(size(tmat,2))));	%indices of nonzero syl-syl transitions for row ii
			% This for loop calculates the numerator of the "loopbackness" score for each trasition type, ranging from 0 to 1
			for k=1:length(thissyl_trans_without)	%indices of nonzero syl-syl transitions for row ii
				loopback1(ii,thissyl_trans_without(k)+2)=starters(thissyl_trans_without(k)+2)*ends(length(ends))+10; %the 10 is just a bias, will be removed
			end
		sylnums=[sylnums,nm(ii)]; %new ev_repeats-style 'nums' vector for only those syllables meeting conditions of this IF loop
	end
end
% loopback1 contains only the numerators of loopback scores for transistions.
% finish loop-corr scores matrix (loopback_scores), create a raw_tmat with loops lumped to endings
loopy_tmat=raw_tmat; total_loopy_trans=0;loopback_scores=zeros(size(raw_tmat,1),size(raw_tmat,2));
findloopback=find(loopback1);	% every transition that exists (except to-begin or to-end), even non-loopbacks
for k=1:length(findloopback)
	loopback_scores(findloopback(k))=(loopback1(findloopback(k))-10)/sum(ends)+10;
	total_loopy_trans=total_loopy_trans+(1-(loopback_scores(findloopback(k))-10));
	[r1,c1]=ind2sub([size(loopback1,1); size(loopback1,2)],findloopback(k));
	loopy_tmat(findloopback(k))=loopy_tmat(findloopback(k))-(loopy_tmat(findloopback(k))*(loopback_scores(findloopback(k))-10));
	loopy_tmat(r1,2)=loopy_tmat(r1,2)+(raw_tmat(findloopback(k))*(loopback_scores(findloopback(k))-10));
end

% Make tmats that are percents
ptmat=tmat.*.01;
for ii=1:size(loopy_tmat,1)
	if sum(loopy_tmat(ii,:))~=0
		loopy_ptmat(ii,:)=loopy_tmat(ii,:)./sum(loopy_tmat(ii,:));
	else
		loopy_ptmat(ii,:)=loopy_tmat(ii,:);
	end
	loopy_ptmat2(ii,:)=loopy_tmat(ii,:)./nm(ii);
end
if loopy_ptmat~=loopy_ptmat2
	fprintf(1,'          ***Incomplete labeling warning\n\r          # of times a syllable is sung\n\r    does not equal sum of observed transitions***\n')  
end

loopy_tmat
raw_tmat

%%%%
% calculations begin here...
total_domtrans=[0,0,0,0];lin_con_sum=[0 0 0 0];entr_sums=[0 0 0 0];
% vector variables: (1)=withendings, nonloopback. (2)=withendings,loopback 3 and 4 are without endings
for ii=1:size(ptmat,1);
	if (~strcmp(nt(ii),'0'))&(~strcmp(nt(ii),'+'))&(~strcmp(nt(ii),'/'));
        	%(~strcmp(nt(ii),'i'))&(~strcmp(nt(ii),'j')
		total_domtrans=[(total_domtrans(1)+max(raw_tmat(ii,:))),(total_domtrans(2)+max(loopy_tmat(ii,:))),(total_domtrans(3)+max(raw_tmat(ii,3:size(raw_tmat,2)))),(total_domtrans(4)+max(loopy_tmat(ii,3:size(loopy_tmat,2))))];	%sums of # dom transitions sung: [withregular,withloopy,withoutregular,withoutloopy]
		% First, non-loopback LC and entr. without endings and with.
		for g = 1:4
			sum_thissyl_lincon=0.0;
			sum_thissyl_entr=0.0;
			if g == 1
			mat=ptmat; s=1; cond='withend, nonloopback'; rmat=raw_tmat;
			end
			if g== 2
			mat=loopy_ptmat; s=1;cond='withend, loopback'; rmat=loopy_tmat;
			end
			if g==3
			mat=ptmat; s=3;cond='exclude ends, nonloopback'; rmat=raw_tmat;
			end
			if g==4
			mat=loopy_ptmat; s=3;cond='exclude ends, loopback'; rmat=loopy_tmat;
			end
			todo=find(mat(ii,s:length(mat(ii,:))));
			if length(todo)<=1
				sum_thissyl_lincon=1;
			% 	so perfectly stereotyped transitions get a '1' for LC and keep the '0' for entropy
			else
				for inrow=1:length(todo);
					sum_thissyl_lincon=sum_thissyl_lincon+(mat(ii,todo(inrow)+(s-1))/(sum(mat(ii,s:length(mat(ii,:))))-mat(ii,todo(inrow)+(s-1))));
					sum_thissyl_entr=sum_thissyl_entr-(mat(ii,todo(inrow)+(s-1))*(log2(mat(ii,todo(inrow)+(s-1)))));
				end
				sum_thissyl_lincon=1-1/sum_thissyl_lincon;
				
			end
			outs(3,g)=outs(3,g) + (sum(rmat(ii,s:size(rmat,2)))/sum(sum(rmat(3:size(rmat,1),s:size(rmat,2)))))*(sum_thissyl_lincon);
			outs(4,g)=outs(4,g)+(sum_thissyl_entr*(sum(rmat(ii,s:length(rmat(ii,:))))/sum(sum(rmat(3:size(rmat,1),s:size(rmat,2))))));
			%outs(3,g)=outs(3,g) + (nm(ii)/sum(sylnums))*(sum_thissyl_lincon);
			%outs(4,g)=outs(4,g)+(sum_thissyl_entr*(nm(ii)/sum(sylnums)));
			if g==condition
				disp_cond=cond; % for display at the end of output
				fprintf(1,['Syllable ''',nt(ii),''': condition = ',cond,'.\n\t\t\btotal transitions = ',num2str(length(find(mat(ii,s:length(mat(ii,:))))))])
				if s==1
			%	if g==1
				if mat(ii,2)~=0
					fprintf(1,[' (*''to-end''*; ',num2str(100*((rmat(ii,2))/sum(rmat(:,2)))),'%% of ends.)'])
					%flags syllables that are song-enders, tells you percent of all endings that this syllable accounts for.
				end
			%	else
			%	if mat(ii,2)~=0
			%		fprintf(1,[' (*''to-end''*; ',num2str(100*nm(ii)*(mat(ii,2)/sum(loopy_tmat(:,2)))),'%% of ends.)'])
			%	end
			%	end
				end
				clear maxval;clear c1;
				[maxval,c1]=max((mat(ii,:)));
				fprintf(1,['\n             Dominant transition = ',num2str(100*maxval),'%%, to ''',num2str(nt(c1)),'.''\n            Transitional entropy = ',num2str(sum_thissyl_entr),'\n              Linear Consistancy = ',num2str(sum_thissyl_lincon),'\n            Prevalence = ',num2str(sum(rmat(ii,s:size(rmat,2)))/sum(sum(rmat(3:size(rmat,1),s:size(rmat,2))))),', Loopback scores = ',mat2str(loopback_scores(ii,find(loopback_scores(ii,:)))-10,3),'\n\n']);
			end
		end
	end
end
% ptmat
fprintf(['  Syllables included for analysis: ',num2str(length(sylnums)),'.\n  Should match total song syls labeled: ',num2str(length(nt)-2),' syllables.\n']);
fprintf(['\n\t\tReg\tLoopbk-corrected','\nDominant\ntrans sung:\t',num2str(total_domtrans(1)),'   \t\t',num2str(total_domtrans(2)),'\nw/o endings:\t',num2str(total_domtrans(3)),'   \t\t',num2str(total_domtrans(4)),',\n\ntotal\ntrans sung:\t',num2str(sum(sum(raw_tmat(3:size(raw_tmat,1),:)))),'   \t\t',num2str(sum(sum(loopy_tmat(3:size(loopy_tmat,1),:)))),'\nw/o endings:\t',num2str(sum(sum(raw_tmat(3:size(raw_tmat,1),3:size(raw_tmat,2))))),'   \t\t',num2str(sum(sum(loopy_tmat(3:size(loopy_tmat,1),3:size(loopy_tmat,2))))),'\n\n'])

%total syl types minus 1, divided by total trans NOT FROM begin or end... = il
outs(1,3)=(size(ptmat,1)-3)/length(find(ptmat(3:size(ptmat,1),3:size(ptmat,2))));
% same as above but use a loopback-corrected transition count... = il loopback-corrected
outs(1,4)=(size(ptmat,1)-3)/total_loopy_trans;
% traditional sl for reference
outs(1,1)=(size(ptmat,1)-2)/length(find(tmat(3:size(tmat,1),:)));
% need to compute a non-internal sl that has loopback correction
outs(1,2) = outs(1,1);
% traditional sc
outs(2,1)=total_domtrans(1)/sum(sum(raw_tmat(3:size(raw_tmat,1),:)));
% sc with loopback correction
outs(2,2)=total_domtrans(2)/sum(sum(loopy_tmat(3:size(loopy_tmat,1),:)));
% sc disregarding to-end transitions-- "iSC"
outs(2,3)=total_domtrans(3)/sum(sum(raw_tmat(3:size(raw_tmat,1),3:size(raw_tmat,2))));
% sc with loopback correction disregarding to-end transitions-- "iSC-loopback-corrected"
outs(2,4)=total_domtrans(4)/sum(sum(loopy_tmat(3:size(loopy_tmat,1),3:size(loopy_tmat,2))));
% JW-style Linear Consistancy, weighted by syllable prevalence. To see unweighted, go to last for loop and swap commented line for real line when defining lin_con_sum
% Unnormalized entropy, weighted by syllable prevalence
%ptmat,loopy_ptmat
%raw_tmat, loopy_tmat
fprintf(1,['\n\t  Regular\tLoopback-Corrected\n\n     SL = ',num2str(outs(1,1)),'\n     IL = ',num2str(outs(1,3)),'      \t\t',num2str(outs(1,4)),'\n\n     SC = ',num2str(outs(2,1)),'      \t\t',num2str(outs(2,2)),'\n    iSC = ',num2str(outs(2,3)),'      \t\t',num2str(outs(2,4)),'\n\n     LC = ',num2str(outs(3,1)),'      \t\t',num2str(outs(3,2)),'\n    iLC = ',num2str(outs(3,3)),'      \t\t',num2str(outs(3,4)),'\n\n   entr = ',num2str(outs(4,1)),'      \t\t',num2str(outs(4,2)),'\n i-entr = ',num2str(outs(4,3)),'      \t\t',num2str(outs(4,4)),'\n'])
%disp(['For the ''',disp_cond,''' condition:'])
lin=outs(1,condition);scon=outs(2,condition);st=(lin+scon)/2;lcon=outs(3,condition);entr=outs(4,condition);
fprintf(1,['\nFor the ''',disp_cond,''' condition:\n','\n       JW_LC = ',num2str(lcon),'\n          ST = ',num2str(st),'\n Entropy sum = ',num2str(entr),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
%disp(['          SL = ',num2str((outs(1,condition)+outs(2,condition))/2)])
%fprintf(1,[Entropy sum = ',num2str(entr_sums(condition)),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
