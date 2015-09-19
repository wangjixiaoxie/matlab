function [sl,sc,st,lc,entr]=newest_test(batch,condition)
%
% cond is an integer 1 - 4. (1)=withendings, nonloopback. (2)=withendings,loopback. 3 and 4 are without endings
% Entropy is reported per-note as total entropy (p*log2(p)), not normalized. Total song-associated entropy
% is reported as sum of each syllable's total entropy (as above) normalized by proportion of times that syllable
% is sung in the batch of songs.
%need:
% pernote: #trans, #non-ending trans--or flag for enders or nonenders, lc value, lc-loopy value, entr, entr norm, entr norm+weighted 
% at end: il, sl, il-loopy, sl-loopy, sc, isc, sc-loopy, lc, lc-loopy, entropy, entr-loopy


[tmat,reps,nt,nm]=ev_repeats(batch,1);

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

%starters= row of tmat for '+' syllable, expressed as percents
starters=[];
starters=.01.*tmat(findstr('+',nt),:);


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
	%	tmat_sum_without=tmat_sum_without+length(thissyl_trans_without);	%this line of code can be moved to main for loop and put next to non-without tmat sum.
		sylnums=[sylnums,nm(ii)]; %new ev_repeats-style 'nums' vector for only those syllables meeting conditions of this IF loop
	end
end
% finish loop-corr scores matrix, create a raw_tmat with loops lumped to endings
loopy_tmat=raw_tmat; total_loopy_trans=0;loopback_scores=zeros(size(raw_tmat,1),size(raw_tmat,2));
findloopback=find(loopback1);	% every transition that exists (except to-begin or to-end), even non-loopbacks
for k=1:length(findloopback)
	loopback_scores(findloopback(k))=(loopback1(findloopback(k))-10)/sum(ends)+10;
	total_loopy_trans=total_loopy_trans+(1-(loopback_scores(findloopback(k))-10));
	[r1,c1]=ind2sub([size(loopback1,1); size(loopback1,2)],findloopback(k));
	loopy_tmat(findloopback(k))=loopy_tmat(findloopback(k))-(loopy_tmat(findloopback(k))*(loopback_scores(findloopback(k))-10));
	loopy_tmat(r1,2)=loopy_tmat(r1,2)+(raw_tmat(findloopback(k))*(loopback_scores(findloopback(k))-10));
end
% tmat, raw_tmat, loopy_tmat

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
		sum_thissyl_lincon=0.0;
		sum_thissyl_entr=0.0;
		% First, non-loopback LC and entr. without endings and with.
		for g = 1:4
			if g == 1
			mat=ptmat; s=1; cond='withend, nonloopback';
			end
			if g== 2
			mat=loopy_ptmat; s=1;cond='withend, loopback';
			end
			if g==3
			mat=ptmat; s=3;cond='exclude ends, nonloopback';
			end
			if g==4
			mat=loopy_ptmat; s=3;cond='exclude ends, loopback';
			end
			todo=find(mat(ii,s:length(mat(ii,:))));
			if length(todo)==1
				sum_thissyl_lincon=1;
			% 	so perfectly stereotyped transitions get a '1' for LC and keep the '0' for entropy
			else
				for inrow=1:length(todo);
					sum_thissyl_lincon=sum_thissyl_lincon+(mat(ii,todo(inrow)+(s-1))/(sum(mat(ii,s:length(mat(ii,:))))-mat(ii,todo(inrow)+(s-1))));
					sum_thissyl_entr=sum_thissyl_entr-(mat(ii,todo(inrow)+(s-1))*(log2(mat(ii,todo(inrow)+(s-1)))));
				end
				sum_thissyl_lincon=1-1/sum_thissyl_lincon;
				
			end
		%	lin_con_sum(1)=lin_con_sum(1) + sum_thissyl_lincon/length(sylnums);
			lin_con_sum(g)=lin_con_sum(g) + (nm(ii)/sum(sylnums))*(sum_thissyl_lincon);
			entr_sums(g)=entr_sums(g)+(sum_thissyl_entr*(nm(ii)/sum(sylnums)));
			if g==condition
				disp_cond=cond; % for display at the end of output
				fprintf(1,['Syllable ''',nt(ii),''': condition = ',cond,'.\n\t\t\btotal transitions = ',num2str(length(find(mat(ii,s:length(mat(ii,:))))))])
				if ptmat(ii,2)~=0
					fprintf(1,[' (** ''to-end'' = ',num2str(100*(raw_tmat(ii,2)/sum(ends))),'%% of ends.)'])
			%		fprintf(1,[' (Including a ''to-end'' transition --\n\t\t\t\t\t\t ''',nt(ii),''' is ',num2str(100*(raw_tmat(ii,2)/sum(ends))),'%% of all song endings.)'])
				end
				clear maxval;clear c1;
				[maxval,c1]=max((tmat(ii,:)));
				fprintf(1,['\n             Dominant transition = ',num2str(maxval),'%%, to ''',num2str(nt(c1)),'.''\n            Transitional entropy = ',num2str(sum_thissyl_entr),'\n              Linear Consistancy = ',num2str(sum_thissyl_lincon),'\n            Prevalence = ',num2str(nm(ii)/sum(sylnums)),', Loopback scores = ',mat2str(loopback_scores(ii,find(loopback_scores(ii,:)))-10,3),'\n\n']);
			end
		end
	end
end
% ptmat
fprintf(['  Syllables included for analysis: ',num2str(length(sylnums)),'.\n  Should match total song syls labeled: ',num2str(length(nt)-2),' syllables.\n']);
fprintf(['\n   Dominant trans sung: ',mat2str(total_domtrans),',\n  and total trans sung: ',mat2str([sum(sum(raw_tmat(3:size(raw_tmat,1),:))),sum(sum(loopy_tmat(3:size(loopy_tmat,1),:))),sum(sum(raw_tmat(3:size(raw_tmat,1),3:size(raw_tmat,2)))),sum(sum(loopy_tmat(3:size(loopy_tmat,1),3:size(loopy_tmat,2))))]),'.\n\n'])

%total syl types minus 1, divided by total trans NOT FROM begin or end
il=(size(ptmat,1)-3)/length(find(ptmat(3:size(ptmat,1),3:size(ptmat,2))));
% same as above but use a loopback-corrected transition count
il_loop=(size(ptmat,1)-3)/total_loopy_trans;
% traditional sl for reference
sl=(size(ptmat,1)-2)/length(find(tmat(3:size(tmat,1),:)));
% traditional sc
sc=total_domtrans(1)/sum(sum(raw_tmat(3:size(raw_tmat,1),:)));
% sc with loopback correction
sc_loop=total_domtrans(2)/sum(sum(loopy_tmat(3:size(loopy_tmat,1),:)));
% sc disregarding to-end transitions-- "iSC"
iSC=total_domtrans(3)/sum(sum(raw_tmat(3:size(raw_tmat,1),3:size(raw_tmat,2))));
% sc with loopback correction disregarding to-end transitions-- "iSC-loopback-corrected"
iSC_loop=total_domtrans(4)/sum(sum(loopy_tmat(3:size(loopy_tmat,1),3:size(loopy_tmat,2))));
% JW-style Linear Consistancy, weighted by syllable prevalence. To see unweighted, go to last for loop and swap commented line for real line when defining lin_con_sum
lc=lin_con_sum(1);lc_loop=lin_con_sum(2);iLC=lin_con_sum(3);iLC_loop=lin_con_sum(4);
% Unnormalized entropy, weighted by syllable prevalence
entr=entr_sums(1);entr_loop=entr_sums(2);int_entr=entr_sums(3);int_entr_loop=entr_sums(4);
fprintf(1,['\tRegular\t\tLoopback-Corrected\n  SL = ',num2str(sl),'\n  IL = ',num2str(il),'\t\t ',num2str(il_loop),'\n  SC = ',num2str(sc),'\t\t ',num2str(sc_loop),'\n iSC = ',num2str(iSC),'\t\t ',num2str(iSC_loop),'\n'])
%disp(['For the ''',disp_cond,''' condition:'])
fprintf(1,['\nFor the ''',disp_cond,''' condition:\n','          SL = ',num2str((sc+sl)/2),'\n       JW_LC = ',num2str(lin_con_sum(condition)),'\n Entropy sum = ',num2str(entr_sums(condition)),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
%fprintf(1,[Entropy sum = ',num2str(entr_sums(condition)),'\n\r\n--------end stereotypy: ',batch,'----------\n\n']);
