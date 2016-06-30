%% lt 9/30/15 - 
% --- specific inputs needed (examples):
% inds=Target_All==0;
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% sign_desired=1; 
% cycles=1000;

%% initiate figure
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];



%% information figure
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

lt_plot_text(0, 0.6, ['desired dir: ' num2str(sign_desired)]);
lt_plot_text(0, 0.5, 'shuffling is done by randomly switching laern and baseline ');

%% collect data

% -- first flip all data to match target learning dir
targ_learn_sign=TargLearnDir_All(inds);
Y_flipped=Y.*repmat(targ_learn_sign',1,2);

% --- Second, match sign of baseline to learning sign (for each syl)
Y_learn_signs=sign(Y_flipped(:,2));
Y_base=Y_flipped(:,1);
Y_base_signmatch=abs(Y_base).*Y_learn_signs;

Y_flipped_BaseSignMatched=[Y_base_signmatch Y_flipped(:,2)];
Y_flipped_BaseSignMatched_original=Y_flipped_BaseSignMatched;

%% Third, take only positive or negative generalizers
if sign_desired==1; % positive
inds_tmp=Y_flipped_BaseSignMatched(:,2)>0;
elseif sign_desired==-1; % negative;
inds_tmp=Y_flipped_BaseSignMatched(:,2)<0;
else
    disp('need a correct sign_desired!! (-1 or 1)');
end
    
% Nx2 vector
Y_flipped_BaseSignMatched=Y_flipped_BaseSignMatched(inds_tmp,:);

% differnces
Y_diffpaired_flipped_BaseSignMatched=Y_flipped_BaseSignMatched(:,2)-Y_flipped_BaseSignMatched(:,1);


%% Perform permutation and overlay plots from all permutation
%% METHOD 1: permuting only using the subset of data that are either negetive or positive generalizers
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('learn minus base(base dir match to learn) shuffled and data overlayed');
xlabel('pitch diff (in targ dir; base sign matched to learning)');
ylabel('cdf');
grid on;

Y_orig=Y_diffpaired_flipped_BaseSignMatched;
N=size(Y_diffpaired_flipped_BaseSignMatched,1);

for i=1:cycles;
% On each round, coin flip for each syl to see if flip baseline and
% learning.
multiplier_shuff=sign(rand(N,1)-0.5);
Y_shuff=multiplier_shuff.*Y_orig;

% plot
[F, X]=ecdf(Y_shuff);
h=plot(X, F, '-', 'Color',[0.5 0.5 0.5],'LineWidth',0.5);
end

% --- PLOT ORIGINAL DATA
lt_plot_cdf(Y_diffpaired_flipped_BaseSignMatched,'b');


% ==== PLOT ORIGINAL DATA ALONE
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('original data');
xlabel('pitch diff (in targ dir; base sign matched to learning)');
ylabel('cdf');
grid on;
lt_plot_cdf(Y_diffpaired_flipped_BaseSignMatched,'b');

% sign rank
p=signrank(Y_flipped_BaseSignMatched(:,2), Y_flipped_BaseSignMatched(:,1), 'tail','right');
lt_plot_pvalue(p, 'one-sided (learn>base) sign rank', 1);

% sign rank
p=signrank(Y_flipped_BaseSignMatched(:,2), Y_flipped_BaseSignMatched(:,1), 'tail','left');
lt_plot_pvalue(p, 'one-sided (learn<base) sign rank', 1);


%% ==== How often does a pitch difference of a given magnitude show up in the shuffled data?
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
 title('How often a pitch difference of given magnitude show up in shuffled data');
% ylabel('cdf of shuffled data');
% xlabel('fract syls with pitch diff in desired dir');
grid on;



%% METHOD 2: permuting as above, but plotting not magnitude, but sign.
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
 title('fract of syls with sign of (learn-base) in desired dir');
ylabel('cdf of shuffled data');
xlabel('fract syls with pitch diff in desired dir');
grid on;


% ORIGINAL DATA
Ysigns=sign(Y_diffpaired_flipped_BaseSignMatched);
Y_InDirOfDesiredSign=sum(Ysigns==sign_desired);
Y_orig=Y_InDirOfDesiredSign/length(Y_diffpaired_flipped_BaseSignMatched);

% SHUFFLED
N=size(Y_diffpaired_flipped_BaseSignMatched,1);
Y_dat_shuff=[];
for i=1:cycles;
% On each round, coin flip for each syl to see if flip baseline and
% learning.
multiplier_shuff=sign(rand(N,1)-0.5);
Y_shuff=multiplier_shuff.*Y_diffpaired_flipped_BaseSignMatched;

% signs
Ysigns_shuff=sign(Y_shuff);
Y_InDirOfDesiredSign_shuff=sum(Ysigns_shuff==sign_desired);

Y_dat_shuff(i)=Y_InDirOfDesiredSign_shuff/length(Y_diffpaired_flipped_BaseSignMatched);
end

lt_plot_cdf(Y_dat_shuff, 'k');
line([Y_orig Y_orig],ylim);








%% METHOD 2: permuting across all the data (i.e. first shuffle, then take the negative or positive ones)
% IN PROGRESS - first shuffle (flp learn and base), then take positive
% "learnign", then do as above.
% ========== 1)first shuffle all data





