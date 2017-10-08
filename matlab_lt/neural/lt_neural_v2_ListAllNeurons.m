function lt_neural_v2_ListAllNeurons
%% LT 2/6/17 - lists all neurons that have been saved to summarystruct

currdir = pwd;
%%
% % load
% SummaryStruct = lt_neural_v2_LoadSummary;
% 
% ==== LIST
[NeuronDatabase] = lt_neural_v2_ConvertSummary2Database(''); % convert to database (then can make table)

% TABLE
NeuronDatabase_table=struct2table(NeuronDatabase.neurons);
summary(NeuronDatabase_table);
disp(NeuronDatabase_table)

positions = get(0, 'MonitorPositions');
f = figure('Position',positions(1,:));
% dat = NeuronDatabase_table;
dat=table2cell(NeuronDatabase_table(:,[1:9 11:end]));
cnames=fieldnames(NeuronDatabase_table(:,[1:9 11:end]));
positions2=get(f, 'Position');
positions2(3:4) = 0.9*positions2(3:4);
rownames = table2cell(NeuronDatabase_table(:,1));
t = uitable('Parent', f, 'Data',dat,'ColumnName',cnames, 'Position', positions2, ...
    'ColumnWidth', {'auto'}, 'RowName', rownames);

% === if there are any empty neurons, then move another neuron

%%
cd(currdir);
