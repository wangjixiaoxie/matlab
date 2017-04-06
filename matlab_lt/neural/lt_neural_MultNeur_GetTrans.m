function [TRANSMATRIX] = lt_neural_MultNeur_GetTrans(NeuronDatabase)

NumNeurons=length(NeuronDatabase.neurons);

maxInterval=0.25; % sec (offset - offset, to count as a transition)



%% 
clear TRANSMATRIX
for i=1:NumNeurons
    try
        cd(NeuronDatabase.global.basedir);
    catch err
        cd(NeuronDatabase.neurons(i).basedir);
    end
        
        % - find day folder
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
    [SongDat, ~, Params] = lt_neural_ExtractDat(batchf, channel_board);

    % ===== EXTRACT TRANSITION MATRIX (i.e. counts)
    % 1) --- get list of all potential syls (including "-");
    allsyls=unique(SongDat.AllLabels);
    allsyls=sort(allsyls); % these are inds for transition matrix    

    % 2) --- initiate transition matrix
    TRANSMATRIX.neuron(i).matrix=zeros(length(allsyls), length(allsyls));
    TRANSMATRIX.neuron(i).sysInOrder=allsyls;
    
    
    % 2) -- fill transition matrix.
    for j=1:length(SongDat.AllLabels)-1
        % slide thru all pairs of syls, based on pre and post syl,
        % increment count in the correct cell in the transition matrix
        inds=[j j+1];
        pairSyls=SongDat.AllLabels(inds);
        offset_to_onset_dur=SongDat.AllOnsets(inds(2))-SongDat.AllOffsets(inds(1));
        
        % --- skip if interval is too long
        if offset_to_onset_dur > maxInterval
            continue
        end
        
        % ---
        ind_firstsyl=strfind(TRANSMATRIX.neuron(i).sysInOrder, pairSyls(1));
        ind_secondsyl=strfind(TRANSMATRIX.neuron(i).sysInOrder, pairSyls(2));
        
        TRANSMATRIX.neuron(i).matrix(ind_firstsyl, ind_secondsyl) = ...
            TRANSMATRIX.neuron(i).matrix(ind_firstsyl, ind_secondsyl) + 1; % increment this cell.
    end

end


%% ==== PLOT AS HEATMAP (WITH NUMBERS)



for i=1:NumNeurons
    
    lt_figure; hold on;
    % --- convert syls to cell
    tmpfun=@(x){x}; 
    cnames=arrayfun(tmpfun, TRANSMATRIX.neuron(i).sysInOrder);

    % ============ frequencies
    lt_subplot(3,1,1); hold on;
    title('frequencies');
    xlabel('to'); ylabel('from');
    dat=TRANSMATRIX.neuron(i).matrix;
    tmp=sort(dat(:)); % ignore - to -
    imagesc(1:size(dat,1), 1:size(dat,2), dat, [0 tmp(end-1)]);
	colormap('autumn');
    set(gca, 'XTick', 1:size(dat,1), 'XTickLabel', cnames, ...
        'YTick', 1:size(dat, 2), 'YTickLabel', cnames);
    % -- put number in each box for each position
    for j=1:size(dat,1);
        for jj=1:size(dat,2);
            dattmp=dat(j, jj); % from j to jj
            lt_plot_text(jj-0.3, j-0.3, num2str(dattmp), 'k');
        end
    end
    % --- get marginals ()
    for j=1:size(dat,1)
        margsum=sum(dat(:,j));
        lt_plot_text(j-0.3, size(dat,1)+0.7, num2str(margsum), 'b');
    end
    % --- get marginals ()
    for j=1:size(dat,1)
        margsum=sum(dat(j,:));
        lt_plot_text(size(dat,1)+0.7, j-0.3, num2str(margsum), 'b');
    end
    axis tight
%     colorbar
    
    % ============ convergent prob
    lt_subplot(3,1,2); hold on;
    title('convergent prob');
    xlabel('to'); ylabel('from');
    
    dat=TRANSMATRIX.neuron(i).matrix;
    dat=dat./repmat(sum(dat, 1), size(dat,1), 1);
    TRANSMATRIX.neuron(i).matrix_convPr=dat;
    
    imagesc(1:size(dat,1), 1:size(dat,2), dat, [0 1]);
%     colorbar
    set(gca, 'XTick', 1:size(dat,1), 'XTickLabel', cnames, ...
        'YTick', 1:size(dat, 2), 'YTickLabel', cnames);
    % -- put number in each box for each position
    for j=1:size(dat,1);
        for jj=1:size(dat,2);
            dattmp=dat(j, jj); % from j to jj
            lt_plot_text(jj-0.3, j-0.3, num2str(100*dattmp, '%2.0f'), 'k');
        end
    end
    axis tight
    
    % ============ divergent prob
    lt_subplot(3,1,3); hold on;
    title('divergent prob');
    xlabel('to'); ylabel('from');
    
    dat=TRANSMATRIX.neuron(i).matrix;
    dat=dat./repmat(sum(dat, 2), 1, size(dat,2));
        TRANSMATRIX.neuron(i).matrix_divPr=dat;

    imagesc(1:size(dat,1), 1:size(dat,2), dat, [0 1]);
%     colorbar
    set(gca, 'XTick', 1:size(dat,1), 'XTickLabel', cnames, ...
        'YTick', 1:size(dat, 2), 'YTickLabel', cnames);
    % -- put number in each box for each position
    for j=1:size(dat,1);
        for jj=1:size(dat,2);
            dattmp=dat(j, jj); % from j to jj
            lt_plot_text(jj-0.3, j-0.3, num2str(100*dattmp, '%2.0f'), 'k');
        end
    end
    axis tight

    
    lt_subtitle(['neuron ' num2str(i)]);
%     
%     
%     
%     
% % TABLES frequencies
% tmp=array2table(TRANSMATRIX.neuron(i).matrix);
% disp(tmp);
% 
% f = figure('Position',[200 200 1000 500]);
% dat=table2cell(tmp);
%     tmpfun=@(x){x}
% cnames=arrayfun(tmpfun,TRANSMATRIX.neuron(i).sysInOrder);
% t = uitable('Parent', f, 'Data',dat,'ColumnName',cnames, 'Position',[0 0 900 450]);
% 
end

%% ========= 


%% ==== troubleshooting
if (0)
% -- plot histogram of offset-onset intervals

tmp=SongDat.AllOnsets(2:end) - SongDat.AllOffsets(1:end-1);
figure; hold on;
lt_plot_histogram(tmp, '', 1, 0, '', 0, 'k');
end

