% Looks at several parameters of wiener entropy per day during learning (mean, sd, and
% cv) using permutation tests.

% Written by DM Brady 10/2012

%% Loads data from bird of interest

which_computer = input('Bluejay number? ');
mult_entropy_boot.parameters.computer = ['/bluejay' num2str(which_computer) '/dbrady/all_wiener_entropy/'];
clear which_computer

mult_entropy_boot.parameters.birdname = input('What is the name of your bird?  ', 's');

mult_entropy_boot.parameters.experiment_type = [input('What type of experiment was it? (pitch, seq, entropy)  ', 's') '/'];

fprintf('\nLoading data...')

try
    load([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type ...
        mult_entropy_boot.parameters.birdname '/' mult_entropy_boot.parameters.birdname '.mat'])
catch err
    display('bird not found, did you run db_wiener_entropy_over_time?')
    return
end

fprintf('done!\n')

%% Loads matlab pool to do parallel processing **note: may take some time**

fprintf('\nSetting up parallel processing...')

%close any existing sessions of matlab
if matlabpool('size') > 0
else
    matlabpool open
end

mult_entropy_boot.parameters.parallel = statset('UseParallel','always');

fprintf('done!\n')

%% Bootstrap mean, sd, cv, median, and iqr over time **note: will take a long time**

fprintf('\nBootstrapping entropy data...')

mult_entropy_boot.parameters.numbootstrap = 10000;
mult_entropy_boot.parameters.names_parameters = {'mean' 'standard deviation' 'coefficient of variation' 'median' 'interquartile range'};
mult_entropy_boot.parameters.units = {'Entropy' 'Entropy' 'no units' 'Entropy' 'Entropy'};

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean)
        try
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j} =...
                bootstrp(mult_entropy_boot.parameters.numbootstrap, @(x) [mean(x) std(x) abs(std(x)./mean(x)*100) median(x) iqr(x)],...
                wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j}(:,2),...
                'Options', mult_entropy_boot.parameters.parallel);
            
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).median{j} =...
                median(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j});
            
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).CI{j} = ...
                [prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j},2.5); ...
                prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j},97.5)];
            
        catch err
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).median{j} ={};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).CI{j} = {};
        end
    end
end

fprintf('done!\n')

%% Calculates baseline

fprintf('\nCalculating baseline...')

% %Calculating baseline (concatenating all entropy before learning into one
% %vector)
% for i = 1:length(wiener_entropy.parameters.syllables)
%     temp_baseline = [];
%     for j = 1:wiener_entropy.parameters.duration.begin
%         temp_baseline = [temp_baseline; wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).time_and_FF{j}(:,2)];
%     end
%     
%     %bootstrapping baseline parameters
%     mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap = ...
%         bootstrp(mult_FF_boot.parameters.numbootstrap, @(x) [mean(x) std(x) std(x)./mean(x) median(x) iqr(x)],...
%         temp_baseline,...
%         'Options', mult_FF_boot.parameters.parallel);
%     
%     mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.median =...
%         median(mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap,1);
%     mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI =...
%         [prctile(mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap,2.5);...
%         prctile(mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap,97.5)];
%     clear temp_baseline
%     
% end

%Calculating baseline (averaging across days)
for i = 1:length(wiener_entropy.parameters.syllables)
    temp_baseline =...
        zeros([mult_entropy_boot.parameters.numbootstrap length(mult_entropy_boot.parameters.names_parameters) wiener_entropy.parameters.duration.begin]);
    for j = 1:wiener_entropy.parameters.duration.begin
        temp_baseline(:,:,j) = mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j};
    end
    
    %calculating baseline parameters
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap = ...
        mean(temp_baseline,3);
    
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.median =...
        median(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap,1);
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI =...
        [prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap,2.5);...
        prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap,97.5)];
    clear temp_baseline
    
end


fprintf('done!\n')   

%% Makes plots of mean, sd, cv, median, and iqr over time

fprintf('\nMaking bootstrap figures...')

if exist([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/bootstrap_figures'], 'dir') ~= 7
    mkdir([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/bootstrap_figures'])
else
end

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(mult_entropy_boot.parameters.names_parameters)
        figure()
        hold on
        grid on
        title([mult_entropy_boot.parameters.names_parameters{j} ' entropy for syllable '...
            wiener_entropy.parameters.syllables(i) ' for ' mult_entropy_boot.parameters.birdname], 'Interpreter', 'None')
        ylabel([mult_entropy_boot.parameters.names_parameters{j} ' (' mult_entropy_boot.parameters.units{j} ')'], 'Interpreter', 'None')
        xlabel('day')
        xlim([-1-wiener_entropy.parameters.duration.begin...
            wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}-wiener_entropy.parameters.duration.begin+1])
        
        %plot baseline
        xlim(get(gca, 'XLim'));
        
        rectangle('Position',[min(get(gca, 'XLim'))+0.2, mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(1,j), ...
            max(get(gca,'XLim'))-min(get(gca, 'XLim'))-0.4, mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(2,j)-...
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(1,j)],...
            'EdgeColor', [.85 .85 .85], 'FaceColor', [.85 .85 .85])
        
%         line(get(gca, 'XLim'),...
%             [mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(1,j) ...
%             mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(1,j)],...
%             'LineStyle', '--', 'Color', [0 0 0])
%         
%         line(get(gca, 'XLim'),...
%             [mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(2,j) ...
%             mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.CI(2,j)],...
%             'LineStyle', '--', 'Color', [0 0 0])
        
        %plots each day's data
        for k = 1:length(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap)
            try
                %plots median from bootstrap
                plot(k-1-wiener_entropy.parameters.duration.begin,...
                     mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).median{k}(j),'o',...
                     'MarkerSize',10,...
                     'MarkerFaceColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k},...
                     'MarkerEdgeColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k})
                 
                 %plots 95% CI from bootstrap
                 line([k-1-wiener_entropy.parameters.duration.begin k-1-wiener_entropy.parameters.duration.begin],...
                     [mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).CI{k}(1,j)...
                     mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).CI{k}(2,j)],...
                     'LineWidth', 2,...
                     'Color', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k})
            catch err
                continue
            end
        end
     
                
        %plots lines for different parts of paradigm
        ylim(get(gca,'YLim'));
        
        try
            line([wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33 ...
                wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([wiener_entropy.parameters.duration.shift-.33 ...
                wiener_entropy.parameters.duration.shift-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33 ...
                wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end

        
        %save figures
        saveas(figure(gcf),[mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
            mult_entropy_boot.parameters.birdname '/bootstrap_figures/' ...
            strrep(mult_entropy_boot.parameters.names_parameters{j},' ','_') '_for_' wiener_entropy.parameters.syllables(i) ...
            '_for_' mult_entropy_boot.parameters.birdname], 'fig')
    end
end

close all

fprintf('done!\n')

%% Calculates % change from baseline

fprintf('\nCalculating percent change from baseline...')

%Calculating percent change from baseline
for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap)
        try
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j} = ...
                (mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j}-...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap)./...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap*100;
            
            %calculating median percent change from baseline
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.median{j} = ...
                median(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j});
            
            %calculating 95% CI percent change from baseline
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.CI{j} =...
                [prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j},2.5);...
                prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j},97.5)];
            
        catch err
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.median{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.CI{j} = {};
            
        end
    end
end

fprintf('done!\n')

%% Makes plot of % change from baseline

fprintf('\nMaking percent change from baseline figures...')

if exist([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/percent_change_figures'], 'dir') ~= 7
    mkdir([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/percent_change_figures'])
else
end

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(mult_entropy_boot.parameters.names_parameters)
        figure()
        hold on
        grid on
        title(['Percent change from baseline: ' mult_entropy_boot.parameters.names_parameters{j} ' for syllable '...
            wiener_entropy.parameters.syllables(i) ' for ' mult_entropy_boot.parameters.birdname], 'Interpreter', 'None')
        ylabel(['Percent change from baseline (%) ' mult_entropy_boot.parameters.names_parameters{j}], 'Interpreter', 'None')
        xlabel('day')
        xlim([-1-wiener_entropy.parameters.duration.begin...
            wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}-wiener_entropy.parameters.duration.begin+1])
        
        %plot baseline
        xlim(get(gca, 'XLim'));
        
        line(get(gca, 'XLim'), [0 0],...
            'LineWidth',2,...
            'LineStyle', '--',...
            'Color', [0 0 0])
        
        %plots each day's data
        for k = 1:length(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap)
            try
                %plots median from bootstrap
                plot(k-1-wiener_entropy.parameters.duration.begin,...
                     mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.median{k}(j),'o',...
                     'MarkerSize',10,...
                     'MarkerFaceColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k},...
                     'MarkerEdgeColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k})
                 
                 %plots 95% CI from bootstrap
                 line([k-1-wiener_entropy.parameters.duration.begin k-1-wiener_entropy.parameters.duration.begin],...
                     [mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.CI{k}(1,j)...
                     mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.CI{k}(2,j)],...
                     'LineWidth', 2,...
                     'Color', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k})
            catch err
                continue
            end
        end
     
                
        %plots lines for different parts of paradigm
        ylim(get(gca,'YLim'));
        
        try
            line([wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33 ...
                wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([wiener_entropy.parameters.duration.shift-.33 ...
                wiener_entropy.parameters.duration.shift-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33 ...
                wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end

        
        %save figures
        saveas(figure(gcf),[mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
            mult_entropy_boot.parameters.birdname '/percent_change_figures/Percent_change_' ...
            strrep(mult_entropy_boot.parameters.names_parameters{j},' ','_') '_for_' wiener_entropy.parameters.syllables(i) ...
            '_for_' mult_entropy_boot.parameters.birdname], 'fig')
    end
end

close all

fprintf('done!\n')

%% Calculates change from baseline

fprintf('\nCalculating change from baseline...')

%Calculating change from baseline
for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap)
        try
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j} = ...
                (mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap{j}-...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.bootstrap);
            
            %calculating median change from baseline
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.median{j} = ...
                median(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j});
            
            %calculating 95% CI change from baseline
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.CI{j} =...
                [prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j},2.5);...
                prctile(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j},97.5)];
            
        catch err
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.median{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.CI{j} = {};
            
        end
    end
end

fprintf('done!\n')


%% Makes plot of change from baseline

fprintf('\nMaking change from baseline figures...')

if exist([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/change_from_baseline_figures'], 'dir') ~= 7
    mkdir([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/change_from_baseline_figures'])
else
end

for i = 1:length(wiener_entropy.parameters.syllables)
    for j = 1:length(mult_entropy_boot.parameters.names_parameters)
        figure()
        hold on
        grid on
        title(['Change from baseline: ' mult_entropy_boot.parameters.names_parameters{j} ' for syllable '...
            wiener_entropy.parameters.syllables(i) ' for ' mult_entropy_boot.parameters.birdname], 'Interpreter', 'None')
        ylabel([mult_entropy_boot.parameters.names_parameters{j} ' (' mult_entropy_boot.parameters.units{j} ')'], 'Interpreter', 'None')
        xlabel('day')
        xlim([-1-wiener_entropy.parameters.duration.begin...
            wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}-wiener_entropy.parameters.duration.begin+1])
        
        %plot baseline
        xlim(get(gca, 'XLim'));
        
        line(get(gca, 'XLim'), [0 0],...
            'LineWidth',2,...
            'LineStyle', '--',...
            'Color', [0 0 0])
        
        %plots each day's data
        for k = 1:length(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).bootstrap)
            try
                %plots median from bootstrap
                plot(k-1-wiener_entropy.parameters.duration.begin,...
                     mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.median{k}(j),'o',...
                     'MarkerSize',10,...
                     'MarkerFaceColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k},...
                     'MarkerEdgeColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k})
                 
                 %plots 95% CI from bootstrap
                 line([k-1-wiener_entropy.parameters.duration.begin k-1-wiener_entropy.parameters.duration.begin],...
                     [mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.CI{k}(1,j)...
                     mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.change_from_baseline.CI{k}(2,j)],...
                     'LineWidth', 2,...
                     'Color', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{k})
            catch err
                continue
            end
        end
     
                
        %plots lines for different parts of paradigm
        ylim(get(gca,'YLim'));
        
        try
            line([wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33 ...
                wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([wiener_entropy.parameters.duration.shift-.33 ...
                wiener_entropy.parameters.duration.shift-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33 ...
                wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end

        
        %save figures
        saveas(figure(gcf),[mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
            mult_entropy_boot.parameters.birdname...
            '/change_from_baseline_figures/change_from_baseline_figures_' strrep(mult_entropy_boot.parameters.names_parameters{j},' ','_')...
            '_for_' wiener_entropy.parameters.syllables(i) '_for_' mult_entropy_boot.parameters.birdname], 'fig')
    end
end

close all

fprintf('done!\n')

%% Making a matrix of percent change from baseline to calculate correlations

fprintf('\nCalculating correlations...')

if exist([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/correlation_figures'], 'dir') ~= 7
    mkdir([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/correlation_figures'])
else
end

for i = 1:length(wiener_entropy.parameters.syllables)
    
    temp_cell_index = find(~cellfun('isempty',mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.median))';
    temp_cell = cell(length(temp_cell_index),1);
    temp_cell_color = cell(length(temp_cell_index),1);
    temp_cell_CI_up = cell(length(temp_cell_index),1);
    temp_cell_CI_lo = cell(length(temp_cell_index),1);
    for j = 1:sum(~cellfun('isempty',mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.median))
        temp_cell{j} = mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.median{temp_cell_index(j)};
        temp_cell_color{j} = wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{temp_cell_index(j)};
        temp_cell_CI_lo{j} = mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.CI{temp_cell_index(j)}(1,:);
        temp_cell_CI_up{j} = mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.percent_from_baseline.CI{temp_cell_index(j)}(2,:);
    end
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat = cell2mat(temp_cell);
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.color = cell2mat(temp_cell_color);
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.CI_lo = cell2mat(temp_cell_CI_lo);
    mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.CI_up = cell2mat(temp_cell_CI_up);
    
    %% calculating correlations between parameters
    [mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.correlations.r...
       mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.correlations.p...
       mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.correlations.rlo...
       mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.correlations.rup] = ...
       corrcoef(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat);
   
   % plot % change in entropy vs. CV
   figure
   hold on
   grid on
   title(['% Change in entropy vs. CV from baseline for syllable ' wiener_entropy.parameters.syllables(i) ' for ' mult_entropy_boot.parameters.birdname])
   xlabel('% Change in entropy from baseline')
   ylabel('% Change in CV from baseline')
%    zlabel('Day')
   
%     plot(mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat(:,1),...
%            mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat(:,3))

%    plot3(mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat(:,1),...
%            mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat(:,3),...
%            temp_cell_index-1-wiener_entropy.parameters.duration.begin)
   
   for j = 1:length(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat)
       
       plot(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,1),...
           mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,3),...
           'o', 'MarkerSize', 10,...
           'MarkerEdgeColor',  mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.color(j,:),...
           'MarkerFaceColor', mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.color(j,:))
       
       line([mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.CI_lo(j,1) ...
           mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.CI_up(j,1)], ...
           [mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,3), ...
           mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,3)],...
           'LineWidth', 2, ...
           'Color', mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.color(j,:))
       
       line([mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,1), ...
           mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,1)],...
           [mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.CI_lo(j,3) ...
           mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.CI_up(j,3)],...
           'LineWidth', 2, ...
           'Color', mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat.color(j,:))
       
%        plot3(mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat(j,1),...
%            mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat(j,3),...
%            temp_cell_index(j)-1-wiener_entropy.parameters.duration.begin,...
%            'o', 'MarkerSize', 10,...
%            'MarkerEdgeColor',  mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat_color(j,:),...
%            'MarkerFaceColor', mult_FF_boot.(wiener_entropy.parameters.syllables(i)).calculations.cell_to_mat_color(j,:))
   end
   
   h.y = get(gca,'Ylim');
   h.x = get(gca,'XLim');
   text(min(h.x)*(3/4),max(h.y)*(3/4),['r = ' ...
       num2str(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.correlations.r(1,3))],...
       'FontSize', 14, 'FontWeight', 'Bold');
   text(min(h.x)*(3/4),max(h.y)*(2/3),['p = ' ...
       num2str(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.correlations.p(1,3))],...
       'FontSize', 14, 'FontWeight', 'Bold');
   
   clear temp*
   clear h
   
   %save figures
   saveas(figure(gcf),[mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
       mult_entropy_boot.parameters.birdname...
       '/correlation_figures/correlation_figures_entropy_vs_CV_for_' ...
       wiener_entropy.parameters.syllables(i) '_for_' mult_entropy_boot.parameters.birdname], 'fig')


   
   
end
 
clear err
clear i
clear j
clear k
close all
    
fprintf('done!\n')  

%% Calculating z score (just for entropy)

fprintf('\nCalculating z scores...')

if exist([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/z_score_figures'], 'dir') ~= 7
    mkdir([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname '/z_score_figures'])
else
end

for i = 1:length(wiener_entropy.parameters.syllables)
    figure
    hold on
    title(['Z score of entropy for syllable '...
        wiener_entropy.parameters.syllables(i) ' for ' mult_entropy_boot.parameters.birdname], 'Interpreter', 'None')
    ylabel(['Z score'], 'Interpreter', 'None')
    xlabel('day')
    xlim([-1-wiener_entropy.parameters.duration.begin...
        wiener_entropy.parameters.days{4}-wiener_entropy.parameters.days{3}-wiener_entropy.parameters.duration.begin+1])
    
    %plot baseline
    xlim(get(gca, 'XLim'));
    
    rectangle('Position',[min(get(gca, 'XLim'))+0.2, -1, ...
        max(get(gca,'XLim'))-min(get(gca, 'XLim'))-0.4, 2],...
        'EdgeColor', [.85 .85 .85], 'FaceColor', [.85 .85 .85])
    
    for j = 1:length(wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean)
        try
            %calculating z score
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.z_score{j} = ...
                (wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).mean{j}(:,2)-...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.median(1))./...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.baseline.median(2);
            
            %calculating mean z score of day
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.mean{j} = ...
                mean(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.z_score{j});
            %calculating sd of z score of day
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.sd{j} = ...
                std(mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.z_score{j});
            
            %plots mean and sd of z score for each day
            plot(j-1-wiener_entropy.parameters.duration.begin,...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.mean{j},'o',...
                'MarkerSize',10,...
                'MarkerFaceColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j},...
                'MarkerEdgeColor', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j})
            
            errorbar(j-1-wiener_entropy.parameters.duration.begin,...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.mean{j},...
                mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.sd{j},...
                'LineWidth', 2,...
                'Color', wiener_entropy.entropy.(wiener_entropy.parameters.syllables(i)).color{j})  
        catch err
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.z_score{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.mean{j} = {};
            mult_entropy_boot.(wiener_entropy.parameters.syllables(i)).calculations.z_score.sd{j} = {};
        end
    end
    
    %plots lines for different parts of paradigm
    ylim(get(gca,'YLim'));
    
    try
        line(get(gca,'XLim'),[0 0],...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    try
        line([wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33 ...
            wiener_entropy.parameters.duration.begin-wiener_entropy.parameters.duration.begin-.33],...
            get(gca, 'YLim'),...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    try
        line([wiener_entropy.parameters.duration.shift-.33 ...
            wiener_entropy.parameters.duration.shift-.33],...
            get(gca, 'YLim'),...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    try
        line([wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33 ...
            wiener_entropy.parameters.duration.shift+wiener_entropy.parameters.duration.consolidation-.33],...
            get(gca, 'YLim'),...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    
    %save figures
    saveas(figure(gcf),[mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
        mult_entropy_boot.parameters.birdname...
        '/z_score_figures/z_score_figure_for_' wiener_entropy.parameters.syllables(i) ...
        '_for_' mult_entropy_boot.parameters.birdname], 'fig')
    
end

close all
clear i
clear j
clear err

fprintf('done!\n') 

%% Save workspace

fprintf('\nSaving workspace...')

save([mult_entropy_boot.parameters.computer mult_entropy_boot.parameters.experiment_type...
    mult_entropy_boot.parameters.birdname...
    '/' mult_entropy_boot.parameters.birdname '_bootstrap.mat'], 'mult_entropy_boot')

    
 fprintf('done!\n')   
    