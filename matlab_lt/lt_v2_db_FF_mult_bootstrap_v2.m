% modified by LT 9/30/13 to fix apaprent error in line 532 where the
% "length" function was giving the size of the incorrect dimension.  

% Looks at several parameters per day during pitch learning (mean, sd, and
% cv) using permutation tests.

% Written by DM Brady 10/2012

%% Loads data from bird of interest

mult_FF_boot.parameters.birdname = input('What is the name of your bird?  ', 's');

mult_FF_boot.parameters.phrase = input('Phrase? ','s');

fprintf('\nLoading data...')

try
    load(['all_days_pitch_' mult_FF_boot.parameters.phrase '/' mult_FF_boot.parameters.birdname '.mat'])
catch err
    display('bird not found, did you run db_contour_and_FF_analysis_over_time?')
    return
end

fprintf('done!\n')

%% Loads matlab pool to do parallel processing **note: may take some time**

% fprintf('\nSetting up parallel processing...')
% 
% % %close any existing sessions of matlab
% % if matlabpool('size') > 0
% % else
% %     matlabpool open
% % end
% 
% %mult_FF_boot.parameters.parallel = statset('UseParallel','always');
% 
% fprintf('done!\n')

%% Bootstrap mean, sd, cv, median, and iqr over time **note: will take a long time**

fprintf('\nBootstrapping FF data...')

mult_FF_boot.parameters.numbootstrap = 10000;
mult_FF_boot.parameters.names_parameters = {'mean' 'standard deviation' 'coefficient of variation' 'median' 'interquartile range'};
mult_FF_boot.parameters.units = {'Hz' 'Hz' 'no units' 'Hz' 'Hz'};

for i = 1:length(multiple_pitch.parameters.syllables)
    for j = 1:length(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF)
        try
%             mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j} =...
%                 bootstrp(mult_FF_boot.parameters.numbootstrap, @(x) [mean(x) std(x) std(x)./mean(x)*100 median(x) iqr(x)],...
%                 multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j}(:,2),...
%                 'Options', mult_FF_boot.parameters.parallel);
            
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j} =...
                bootstrp(mult_FF_boot.parameters.numbootstrap, @(x) [mean(x) std(x) std(x)./mean(x)*100 median(x) iqr(x)],...
                multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j}(:,2));
            
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).median{j} =...
                median(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j});
            
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).CI{j} = ...
                [prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j},2.5); ...
                prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j},97.5)];
            
        catch err
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).median{j} ={};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).CI{j} = {};
        end
    end
end

fprintf('done!\n')

%% Calculates baseline

fprintf('\nCalculating baseline...')

% %Calculating baseline (concatenating all FF before learning into one
% %vector)
% for i = 1:length(multiple_pitch.parameters.syllables)
%     temp_baseline = [];
%     for j = 1:multiple_pitch.parameters.duration.begin
%         temp_baseline = [temp_baseline; multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j}(:,2)];
%     end
%     
%     %bootstrapping baseline parameters
%     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap = ...
%         bootstrp(mult_FF_boot.parameters.numbootstrap, @(x) [mean(x) std(x) std(x)./mean(x) median(x) iqr(x)],...
%         temp_baseline,...
%         'Options', mult_FF_boot.parameters.parallel);
%     
%     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.median =...
%         median(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap,1);
%     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI =...
%         [prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap,2.5);...
%         prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap,97.5)];
%     clear temp_baseline
%     
% end

%Calculating baseline (averaging across days)
for i = 1:length(multiple_pitch.parameters.syllables)
    temp_baseline =...
        zeros([mult_FF_boot.parameters.numbootstrap length(mult_FF_boot.parameters.names_parameters) multiple_pitch.parameters.duration.begin]);
    for j = 1:multiple_pitch.parameters.duration.begin
        temp_baseline(:,:,j) = mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j};
    end
    
    %calculating baseline parameters
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap = ...
        mean(temp_baseline,3);
    
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.median =...
        median(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap,1);
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI =...
        [prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap,2.5);...
        prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap,97.5)];
    clear temp_baseline
    
end


fprintf('done!\n')   

%% Makes plots of mean, sd, cv, median, and iqr over time

fprintf('\nMaking bootstrap figures...')

if exist(['all_days_pitch_' mult_FF_boot.parameters.phrase '/bootstrap_figures'], 'dir') ~= 7
    mkdir(['all_days_pitch_' mult_FF_boot.parameters.phrase '/bootstrap_figures'])
else
end

for i = 1:length(multiple_pitch.parameters.syllables)
    for j = 1:length(mult_FF_boot.parameters.names_parameters)
        figure()
        hold on
        grid on
        title([mult_FF_boot.parameters.names_parameters{j} ' FF for syllable '...
            multiple_pitch.parameters.syllables(i) ' for ' mult_FF_boot.parameters.birdname], 'Interpreter', 'None')
        ylabel([mult_FF_boot.parameters.names_parameters{j} ' (' mult_FF_boot.parameters.units{j} ')'], 'Interpreter', 'None')
        xlabel('day')
        xlim([-1-multiple_pitch.parameters.duration.begin...
            multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}-multiple_pitch.parameters.duration.begin+1])
        
        %plot baseline
        xlim(get(gca, 'XLim'));
        
        rectangle('Position',[min(get(gca, 'XLim'))+0.2, mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(1,j), ...
            max(get(gca,'XLim'))-min(get(gca, 'XLim'))-0.4, mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(2,j)-...
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(1,j)],...
            'EdgeColor', [.85 .85 .85], 'FaceColor', [.85 .85 .85])
        
%         line(get(gca, 'XLim'),...
%             [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(1,j) ...
%             mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(1,j)],...
%             'LineStyle', '--', 'Color', [0 0 0])
%         
%         line(get(gca, 'XLim'),...
%             [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(2,j) ...
%             mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.CI(2,j)],...
%             'LineStyle', '--', 'Color', [0 0 0])
        
        %plots each day's data
        for k = 1:length(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap)
            try
                %plots median from bootstrap
                plot(k-1-multiple_pitch.parameters.duration.begin,...
                     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).median{k}(j),'o',...
                     'MarkerSize',10,...
                     'MarkerFaceColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k},...
                     'MarkerEdgeColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k})
                 
                 %plots 95% CI from bootstrap
                 line([k-1-multiple_pitch.parameters.duration.begin k-1-multiple_pitch.parameters.duration.begin],...
                     [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).CI{k}(1,j)...
                     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).CI{k}(2,j)],...
                     'LineWidth', 2,...
                     'Color', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k})
            catch err
                continue
            end
        end
     
                
        %plots lines for different parts of paradigm
        ylim(get(gca,'YLim'));
        
        try
            line([multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33 ...
                multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([multiple_pitch.parameters.duration.shift-.33 ...
                multiple_pitch.parameters.duration.shift-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33 ...
                multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end

        
        %save figures
        saveas(figure(gcf),['all_days_pitch_' mult_FF_boot.parameters.phrase '/bootstrap_figures/' ...
            strrep(mult_FF_boot.parameters.names_parameters{j},' ','_') '_for_' multiple_pitch.parameters.syllables(i) ...
            '_for_' mult_FF_boot.parameters.birdname], 'fig')
    end
end

close all

fprintf('done!\n')

%% Calculates % change from baseline

fprintf('\nCalculating percent change from baseline...')

%Calculating percent change from baseline
for i = 1:length(multiple_pitch.parameters.syllables)
    for j = 1:length(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap)
        try
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j} = ...
                (mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j}-...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap)./...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap*100;
            
            %calculating median percent change from baseline
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.median{j} = ...
                median(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j});
            
            %calculating 95% CI percent change from baseline
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.CI{j} =...
                [prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j},2.5);...
                prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j},97.5)];
            
        catch err
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.bootstrap{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.median{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.CI{j} = {};
            
        end
    end
end

fprintf('done!\n')

%% Makes plot of % change from baseline

fprintf('\nMaking percent change from baseline figures...')

if exist(['all_days_pitch_' mult_FF_boot.parameters.phrase '/percent_change_figures'], 'dir') ~= 7
    mkdir(['all_days_pitch_' mult_FF_boot.parameters.phrase '/percent_change_figures'])
else
end

for i = 1:length(multiple_pitch.parameters.syllables)
    for j = 1:length(mult_FF_boot.parameters.names_parameters)
        figure()
        hold on
        grid on
        title(['Percent change from baseline: ' mult_FF_boot.parameters.names_parameters{j} ' for syllable '...
            multiple_pitch.parameters.syllables(i) ' for ' mult_FF_boot.parameters.birdname], 'Interpreter', 'None')
        ylabel(['Percent change from baseline (%) ' mult_FF_boot.parameters.names_parameters{j}], 'Interpreter', 'None')
        xlabel('day')
        xlim([-1-multiple_pitch.parameters.duration.begin...
            multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}-multiple_pitch.parameters.duration.begin+1])
        
        %plot baseline
        xlim(get(gca, 'XLim'));
        
        line(get(gca, 'XLim'), [0 0],...
            'LineWidth',2,...
            'LineStyle', '--',...
            'Color', [0 0 0])
        
        %plots each day's data
        for k = 1:length(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap)
            try
                %plots median from bootstrap
                plot(k-1-multiple_pitch.parameters.duration.begin,...
                     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.median{k}(j),'o',...
                     'MarkerSize',10,...
                     'MarkerFaceColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k},...
                     'MarkerEdgeColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k})
                 
                 %plots 95% CI from bootstrap
                 line([k-1-multiple_pitch.parameters.duration.begin k-1-multiple_pitch.parameters.duration.begin],...
                     [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.CI{k}(1,j)...
                     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.CI{k}(2,j)],...
                     'LineWidth', 2,...
                     'Color', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k})
            catch err
                continue
            end
        end
     
                
        %plots lines for different parts of paradigm
        ylim(get(gca,'YLim'));
        
        try
            line([multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33 ...
                multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([multiple_pitch.parameters.duration.shift-.33 ...
                multiple_pitch.parameters.duration.shift-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33 ...
                multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end

        
        %save figures
        saveas(figure(gcf),['all_days_pitch_' mult_FF_boot.parameters.phrase '/percent_change_figures/Percent_change_' ...
            strrep(mult_FF_boot.parameters.names_parameters{j},' ','_') '_for_' multiple_pitch.parameters.syllables(i) ...
            '_for_' mult_FF_boot.parameters.birdname], 'fig')
    end
end

close all

fprintf('done!\n')

%% Calculates change from baseline

fprintf('\nCalculating change from baseline...')

%Calculating change from baseline
for i = 1:length(multiple_pitch.parameters.syllables)
    for j = 1:length(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap)
        try
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j} = ...
                (mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap{j}-...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.bootstrap);
            
            %calculating median change from baseline
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.median{j} = ...
                median(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j});
            
            %calculating 95% CI change from baseline
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.CI{j} =...
                [prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j},2.5);...
                prctile(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j},97.5)];
            
        catch err
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.bootstrap{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.median{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.CI{j} = {};
            
        end
    end
end

fprintf('done!\n')


%% Makes plot of change from baseline

fprintf('\nMaking change from baseline figures...')

if exist(['all_days_pitch_' mult_FF_boot.parameters.phrase '/change_from_baseline_figures'], 'dir') ~= 7
    mkdir(['all_days_pitch_' mult_FF_boot.parameters.phrase '/change_from_baseline_figures'])
else
end

for i = 1:length(multiple_pitch.parameters.syllables)
    for j = 1:length(mult_FF_boot.parameters.names_parameters)
        figure()
        hold on
        grid on
        title(['Change from baseline: ' mult_FF_boot.parameters.names_parameters{j} ' for syllable '...
            multiple_pitch.parameters.syllables(i) ' for ' mult_FF_boot.parameters.birdname], 'Interpreter', 'None')
        ylabel([mult_FF_boot.parameters.names_parameters{j} ' (' mult_FF_boot.parameters.units{j} ')'], 'Interpreter', 'None')
        xlabel('day')
        xlim([-1-multiple_pitch.parameters.duration.begin...
            multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}-multiple_pitch.parameters.duration.begin+1])
        
        %plot baseline
        xlim(get(gca, 'XLim'));
        
        line(get(gca, 'XLim'), [0 0],...
            'LineWidth',2,...
            'LineStyle', '--',...
            'Color', [0 0 0])
        
        %plots each day's data
        for k = 1:length(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).bootstrap)
            try
                %plots median from bootstrap
                plot(k-1-multiple_pitch.parameters.duration.begin,...
                     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.median{k}(j),'o',...
                     'MarkerSize',10,...
                     'MarkerFaceColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k},...
                     'MarkerEdgeColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k})
                 
                 %plots 95% CI from bootstrap
                 line([k-1-multiple_pitch.parameters.duration.begin k-1-multiple_pitch.parameters.duration.begin],...
                     [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.CI{k}(1,j)...
                     mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.change_from_baseline.CI{k}(2,j)],...
                     'LineWidth', 2,...
                     'Color', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{k})
            catch err
                continue
            end
        end
     
                
        %plots lines for different parts of paradigm
        ylim(get(gca,'YLim'));
        
        try
            line([multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33 ...
                multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([multiple_pitch.parameters.duration.shift-.33 ...
                multiple_pitch.parameters.duration.shift-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end
        
        try
            line([multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33 ...
                multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33],...
                get(gca, 'YLim'),...
                'LineStyle', '--', 'Color', [0 0 0])
        catch err
            continue
        end

        
        %save figures
        saveas(figure(gcf),['all_days_pitch_' mult_FF_boot.parameters.phrase ...
            '/change_from_baseline_figures/change_from_baseline_figures_' strrep(mult_FF_boot.parameters.names_parameters{j},' ','_')...
            '_for_' multiple_pitch.parameters.syllables(i) '_for_' mult_FF_boot.parameters.birdname], 'fig')
    end
end

close all

fprintf('done!\n')

%% Making a matrix of percent change from baseline to calculate correlations

fprintf('\nCalculating correlations...')

if exist(['all_days_pitch_' mult_FF_boot.parameters.phrase  '/correlation_figures'], 'dir') ~= 7
    mkdir(['all_days_pitch_' mult_FF_boot.parameters.phrase '/correlation_figures'])
else
end

for i = 1:length(multiple_pitch.parameters.syllables)
    
    temp_cell_index = find(~cellfun('isempty',mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.median))';
    temp_cell = cell(length(temp_cell_index),1);
    temp_cell_color = cell(length(temp_cell_index),1);
    temp_cell_CI_up = cell(length(temp_cell_index),1);
    temp_cell_CI_lo = cell(length(temp_cell_index),1);
    for j = 1:sum(~cellfun('isempty',mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.median))
        temp_cell{j} = mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.median{temp_cell_index(j)};
        temp_cell_color{j} = multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{temp_cell_index(j)};
        temp_cell_CI_lo{j} = mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.CI{temp_cell_index(j)}(1,:);
        temp_cell_CI_up{j} = mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.percent_from_baseline.CI{temp_cell_index(j)}(2,:);
    end
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat = cell2mat(temp_cell);
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.color = cell2mat(temp_cell_color);
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.CI_lo = cell2mat(temp_cell_CI_lo);
    mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.CI_up = cell2mat(temp_cell_CI_up);
    
    %% calculating correlations between parameters
    [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.correlations.r...
       mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.correlations.p...
       mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.correlations.rlo...
       mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.correlations.rup] = ...
       corrcoef(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat);
   
   % plot % change in FF vs. CV
   figure
   hold on
   grid on
   title(['% Change in FF vs. CV from baseline for syllable ' multiple_pitch.parameters.syllables(i) ' for ' mult_FF_boot.parameters.birdname])
   xlabel('% Change in FF from baseline')
   ylabel('% Change in CV from baseline')
%    zlabel('Day')
   
%     plot(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat(:,1),...
%            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat(:,3))

%    plot3(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat(:,1),...
%            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat(:,3),...
%            temp_cell_index-1-multiple_pitch.parameters.duration.begin)
   
   for j = 1:size(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat,1)
       
       plot(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,1),...
           mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,3),...
           'o', 'MarkerSize', 10,...
           'MarkerEdgeColor',  mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.color(j,:),...
           'MarkerFaceColor', mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.color(j,:))
       
       line([mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.CI_lo(j,1) ...
           mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.CI_up(j,1)], ...
           [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,3), ...
           mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,3)],...
           'LineWidth', 2, ...
           'Color', mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.color(j,:))
       
       line([mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,1), ...
           mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.cell_to_mat(j,1)],...
           [mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.CI_lo(j,3) ...
           mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.CI_up(j,3)],...
           'LineWidth', 2, ...
           'Color', mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat.color(j,:))
       
%        plot3(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat(j,1),...
%            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat(j,3),...
%            temp_cell_index(j)-1-multiple_pitch.parameters.duration.begin,...
%            'o', 'MarkerSize', 10,...
%            'MarkerEdgeColor',  mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat_color(j,:),...
%            'MarkerFaceColor', mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.cell_to_mat_color(j,:))
   end
   
   h.y = get(gca,'Ylim');
   h.x = get(gca,'XLim');
   text(min(h.x)*(3/4),max(h.y)*(3/4),['r = ' ...
       num2str(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.correlations.r(1,3))],...
       'FontSize', 14, 'FontWeight', 'Bold');
   text(min(h.x)*(3/4),max(h.y)*(2/3),['p = ' ...
       num2str(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.correlations.p(1,3))],...
       'FontSize', 14, 'FontWeight', 'Bold');
   
   clear temp*
   clear h
   
   %save figures
   saveas(figure(gcf),['all_days_pitch_' mult_FF_boot.parameters.phrase ...
       '/correlation_figures/correlation_figures_FF_vs_CV_for_' ...
       multiple_pitch.parameters.syllables(i) '_for_' mult_FF_boot.parameters.birdname], 'fig')


   
   
end
 
clear err
clear i
clear j
clear k
close all
    
fprintf('done!\n')  

%% Calculating z score (just for FF)

fprintf('\nCalculating z scores...')

if exist(['all_days_pitch_' mult_FF_boot.parameters.phrase '/z_score_figures'], 'dir') ~= 7
    mkdir(['all_days_pitch_' mult_FF_boot.parameters.phrase '/z_score_figures'])
else
end

for i = 1:length(multiple_pitch.parameters.syllables)
    figure
    hold on
    title(['Z score of FF for syllable '...
        multiple_pitch.parameters.syllables(i) ' for ' mult_FF_boot.parameters.birdname], 'Interpreter', 'None')
    ylabel(['Z score'], 'Interpreter', 'None')
    xlabel('day')
    xlim([-1-multiple_pitch.parameters.duration.begin...
        multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}-multiple_pitch.parameters.duration.begin+1])
    
    %plot baseline
    xlim(get(gca, 'XLim'));
    
    rectangle('Position',[min(get(gca, 'XLim'))+0.2, -1, ...
        max(get(gca,'XLim'))-min(get(gca, 'XLim'))-0.4, 2],...
        'EdgeColor', [.85 .85 .85], 'FaceColor', [.85 .85 .85])
    
    for j = 1:length(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF)
        try
            %calculating z score
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.z_score{j} = ...
                (multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j}(:,2)-...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.median(1))./...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.baseline.median(2);
            
            %calculating mean z score of day
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.mean{j} = ...
                mean(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.z_score{j});
            %calculating sd of z score of day
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.sd{j} = ...
                std(mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.z_score{j});
            
            %plots mean and sd of z score for each day
            plot(j-1-multiple_pitch.parameters.duration.begin,...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.mean{j},'o',...
                'MarkerSize',10,...
                'MarkerFaceColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j},...
                'MarkerEdgeColor', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j})
            
            errorbar(j-1-multiple_pitch.parameters.duration.begin,...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.mean{j},...
                mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.sd{j},...
                'LineWidth', 2,...
                'Color', multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).color{j})  
        catch err
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.z_score{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.mean{j} = {};
            mult_FF_boot.(multiple_pitch.parameters.syllables(i)).calculations.z_score.sd{j} = {};
        end
    end
    
    %plots lines for different parts of paradigm
    ylim(get(gca,'YLim'));
    
    try
        line(get(gca,'XLim'), [0 0],...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    try
        line([multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33 ...
            multiple_pitch.parameters.duration.begin-multiple_pitch.parameters.duration.begin-.33],...
            get(gca, 'YLim'),...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    try
        line([multiple_pitch.parameters.duration.shift-.33 ...
            multiple_pitch.parameters.duration.shift-.33],...
            get(gca, 'YLim'),...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    try
        line([multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33 ...
            multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.33],...
            get(gca, 'YLim'),...
            'LineStyle', '--', 'Color', [0 0 0])
    catch err
        continue
    end
    
    
    %save figures
    saveas(figure(gcf),['all_days_pitch_' mult_FF_boot.parameters.phrase ...
        '/z_score_figures/z_score_figure_for_' multiple_pitch.parameters.syllables(i) ...
        '_for_' mult_FF_boot.parameters.birdname], 'fig')
    
end

close all
clear i
clear j
clear err

fprintf('done!\n') 

%% Save workspace

fprintf('\nSaving workspace...')

save(['all_days_pitch_' mult_FF_boot.parameters.phrase ...
    '/' mult_FF_boot.parameters.birdname '_bootstrap.mat'], 'mult_FF_boot')

    
 fprintf('done!\n')   
    
    
    
    
    
    
    
    
    
    






