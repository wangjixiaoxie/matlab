function [lin,scon,st,lcon,entr]=jw_stereot2(batch,condition,INTRO)
% [lin,scon,st,lcon,entr]=jw_stereot2(batch,condition,intro)
%
% condition = integer 1 thru 4.	(1)=withendings, nonloopback.
%				(2)=withendings,loopback-corrected.
%				3 and 4 are without endings
% 
% intro is a syllable label that represents an intronote, for the
% purpose of determining loopback transitions.
%   (default is the + label, which starts each file. Loopback score
%   for a transition depends both on the "ending propensity" for the
%   "from" syllable and on the "starting propensity" of the "to"
%   syllable of the transition in question. So "starting propensity"
%   can be judged based on either how often the song begins with a
%   certain label or how often the POST-INTRONOTE song begins with
%   that label. For loopback correction, songs are considered to
%   start with whatever label follows the "intro" argument.)
%
% Many stereotypy values are computed and reported in standard output,
% but only variables that fit the condition argument returned.
% 
% lin = published 'sequence linearity' (with endings) or internal
%	linearity (excludes endings) according to 'condition' arg.
%	see Scharff & Nottebohm 1991, Foster & Bottjer 2000
% scon= published 'sequence consistancy.'  see Scharff & Nottebohm 1991.
% st  = average of lin and scon, see Scharff & Nottebohm 1991.
% lcon= JW's 'linear consistancy' measure. scores each syllable as
%	having a linear consistancy of 0 - 1. If the transitions
%	for a syl are .99 and .01 probabilities, lc is almost 1,
%	and increasing either the prevalence of the minor trans or
%	the # possible trans will smoothly decrease the score.
%	lc = 1 - 1/(sum(p/(1-p))).  Weighted by prevalence of syls.
% entr= entropy. Weighted sum of entropy (p*log2(p)) for each syllable.
%	weighted by prevalence of syls.
%
% At end of function's output, values reported with or without 'to-end'
% transitions and with or without loopback-correction. 
%
% Per-syllable values are reported first though, wherein only the values
% pertaining to the 'condition' argument are displayed.
% LC is reported PER-NOTE as LC value, not normalized. TOTAL song LC
% (at end of output) is each syl's value, prevalence normalized.
% Entropy is reported PER-NOTE as entropy (p*log2(p)), not
% normalized. TOTAL song-associated entropy is prevalence normalized.
%
% Requires ev_repeats.m.
if (~exist('INTRO'))
	INTRO='+';
end
[tmat,reps,nt,nm]=ev_repeats(batch,1);

outs=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];	%rows for linearity, consistancy, lc, and entr.
					%columns = the 4 conditions

% b= # files listed in batch file
fid=fopen(batch,'r');
b=0;
while (1);
    ln = fgetl(fid);
    if (~ischar(ln));
        break;
    end
    b=b+1;
end
fprintf(1,['\n----------Stereotypy: ',batch,'------------\n\r',num2str(b),' songs in batchfile.\n\r']);fclose(fid);
if findstr(INTRO,nt) ==[]
	fprintf(1,['The requested intronote ''',nt(findstr(INTRO,nt)),''' does not occur.\nUsing the true starts of songs as loop-starters.\n'])
	INTRO='+';
end

%starters= row of tmat for 'intro' syllable, expressed as percents, with 100% being all NON-i-TO-i transitions only
starters=[];
starters=.01.*tmat(findstr(INTRO,nt),:).*nm(findstr(INTRO,nt));	%'raw starters'
starters(findstr(INTRO,nt))=0;
starters=starters./sum(starters);

% Begin to make loopback-correction scores matrix tmat, count good syllables..
ends=[];loopback1=[];sylnums=[];
for ii=1:size(tmat,1);
	% a non-percentage-based tmat:
	raw_tmat(ii,:)=.01.*tmat(ii,:).*nm(ii);
	if (~strcmp(nt(ii),'0'))&(~strcmp(nt(ii),'+'))&(~strcmp(nt(ii),'/'));
			%(~strcmp(nt(ii),'i'))&(~strcmp(nt(ii),'j')
		ends=[ends,(raw_tmat(ii,2))]; %how many song endings does each syllable do? sum is total # of endings..
		thissyl_trans_without=find(tmat(ii,3:(size(tmat,2))));	%indices of nonzero syl-syl transitions for row ii
			% This for loop calculates the numerator of the "loopbackness" score for each trasition type, ranging from 0 to 1
			for k=1:length(thissyl_trans_without)	%indices of nonzero syl-syl transitions for row ii
				loopback1(ii,thissyl_trans_without(k)+2)=starters(thissyl_trans_without(k)+2)*ends(length(ends))+10;
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
