function lt_neural_MultNeur_SpkWaveforms(NeuronDatabase, numrandspikes)

%% ==== plot  spike waveforms for all neurons
close all;
% numrandspikes=100;
NumNeurons=length(NeuronDatabase.neurons);

figcount=1;
subplotrows=NumNeurons;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

plotcols=lt_make_plot_colors(NumNeurons, 0, 0);

for i=1:NumNeurons
    
    cd(NeuronDatabase.neurons(i).basedir)

    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    if length(tmp)>1
        tmp = dir([dirdate '_' NeuronDatabase.neurons(i).exptID '*']);
        assert(length(tmp) ==1,' daiosfhasiohfioawe');
    end
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    
    inds=find(NeurDat.spikes_cat.cluster_class(:,1)==NeuronDatabase.neurons(i).clustnum); % get desired clust
    % get random subsamp
    inds=inds(randperm(length(inds), numrandspikes)); % get subset of spikes
    spkwaves=NeurDat.spikes_cat.spikes(inds, :);
    fs=NeurDat.metaDat(1).fs;
    
    % -- plot individual spikes
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['neuron ' num2str(i)]);
    tt=1:size(spkwaves,2);
    tt=1000*tt/fs;
    plot(tt, spkwaves', 'Color', plotcols{i});
    
    % -- overlay mean + std
    spkwaves_mean=mean(spkwaves,1);
    spkwaves_std=std(spkwaves,0, 1);
    plot(tt, spkwaves_mean, 'Color', plotcols{i}, 'LineWidth', 3);
    
    if strcmp(NeuronDatabase.neurons(i).NOTE_is_single_unit, 'no')
        % then is mu
        lt_plot_text(tt(1), 200, 'MU');
        
    end
    
    axis tight
end
