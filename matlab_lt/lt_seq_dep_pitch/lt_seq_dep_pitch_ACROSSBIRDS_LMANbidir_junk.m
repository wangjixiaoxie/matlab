
%% % ==== SUBPLOT 1 - Before start 2 dir learning [AFP BIAS]
count=1;
SubplotsPerFig=8;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP BIAS [During single WN targ]')
ylabel('pitch shift, norm to target learning (PBS)');

plotcols=lt_make_plot_colors(ExptCounter-1, 0, 0);
Yexptnums={Y_ExptNum_Target, Y_ExptNum_NewTarg, Y_ExptNum_AllOthers};

ind=1; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymeans=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_AFP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_AFP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- All nontargets
x=3;
DatArray=Y_AFP_AllOthers;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
X=1:3;
hbar=lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




%% ==== SUBPLOT 2 - DURING BIDIR (late (early if no late)) [afp bias]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP BIAS [During bidir WN targ]')
ylabel('pitch shift, norm to target learning (PBS)');

ind=2; % during bidir

Yraw={}; % each cell is one category of syllable
Ymeans=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_AFP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_AFP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- All nontargets
x=3;
DatArray=Y_AFP_AllOthers;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
% for i=1:length(Yraw);
%     plot(i+0.2, Yraw{i}, 'ok');
% end

for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end


% 2) means
X=1:3;
hbar=lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);





%% ==== SUBPLOT 3 - Before start 2 dir learning [LEARNING]
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('LEARNING [During single WN targ]')
ylabel('pitch shift, norm to target learning (PBS)');

ind=1; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymeans=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_Learning_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_Learning_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- All nontargets
x=3;
DatArray=Y_Learning_AllOthers;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
X=1:3;
hbar=lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




%% ==== SUBPLOT 4 -
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('LEARNING [During bidir WN targ]')
ylabel('pitch shift, norm to target learning (PBS)');

ind=2; % during bidir

Yraw={}; % each cell is one category of syllable
Ymeans=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_Learning_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_Learning_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- All nontargets
x=3;
DatArray=Y_Learning_AllOthers;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
X=1:3;
hbar=lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);





%% ==== SUBPLOT 5
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP BIAS [During single WN targ]')
ylabel('pitch shift, norm to target learning (PBS)');

ind=1; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymeans=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_MP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_MP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- All nontargets
x=3;
DatArray=Y_MP_AllOthers;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
% 1) raw vals
for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
X=1:3;
hbar=lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);




%% ==== SUBPLOT 5
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('MP BIAS [During bidir WN targ]')
ylabel('pitch shift, norm to target learning (PBS)');

ind=2; % pre-start BIDIR

Yraw={}; % each cell is one category of syllable
Ymeans=[];
Ysem=[];

% --- First targets
x=1;
DatArray=Y_MP_Target;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- Second targets
x=2;
DatArray=Y_MP_NewTarg;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));

% --- All nontargets
x=3;
DatArray=Y_MP_AllOthers;

Yraw{x}=DatArray(:,ind);
Ymean(x)=mean(DatArray(:,ind));
Ysem(x)=lt_sem(DatArray(:,ind));


% ------------------ PLOT BARS
for i=1:length(Yraw);
    for ii=1:ExptCounter-1;
    indsinds=Yexptnums{i}==ii;
    plot(i+0.2, Yraw{i}(indsinds), 'o' ,'Color',plotcols{ii});
    end
end

% 2) means
X=1:3;
hbar=lt_plot_bar(X, Ymean, {'Errors', Ysem, 'Color' 'k'});

Xlabel={'First Targ', 'Second Targ', 'Nontargets (others)'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);

rotateXLabels(gca, 45);
