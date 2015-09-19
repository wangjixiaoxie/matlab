function [all_days_all_analysis]= lt_all_days_PLOT_AllTrans(baseline_days, all_trans_index, con_or_div,motif);
%% LT 3/28/14
% plots transition matrices across days
% run this function from bird folder. must first run lt_all_days_various and then lt_all_days_all_analysis to compile data

% INPUT PARAMETERS(and examples)
% baseline_days = [1 2 3] for 3 days baseline (e.g. for WN learning)
% all_trans_index= "28Mar2014_1925" (date and time of the structure to load)
% con_or_div = "convergence" or "divergence"

close all

% load data structure
cd all_days_all_analysis/
fname=dir(['all_days_all_analysis_' all_trans_index '*']);
load(fname.name)
cd ..

% Figure out the syllables that were sang every day.
days_with_analysis=[];
for i=1:size(all_days_all_analysis.data,2);
    try
        if ~isempty(all_days_all_analysis.data{i}.all_songs.all_trans);
            days_with_analysis=[days_with_analysis i];
        end
    catch err
        continue
    end
end


% concatenate syl labels and data matrices
num_of_days=length(days_with_analysis);
for j=1:num_of_days;
    syl_labels_concatenated{j}=all_days_all_analysis.data{days_with_analysis(j)}.all_songs.all_trans.summary.syl_labels_in_order;
end

% % first get indices for transitions that were present on baseline days.
% clear dummy
% for j=1:length(baseline_days);
%     dummy(:,:,j)=all_days_all_analysis.data{baseline_days(j)}.all_songs.all_trans.summary.(con_or_div).matrix_of_fractions(syl_indices{j},syl_indices{j});
% end
% all_days_all_analysis.PLOT.all_trans.(con_or_div).baseline=mean(dummy,3);
% all_days_all_analysis.PLOT.all_trans.(con_or_div).baseline_extant_transitions=(all_days_all_analysis.PLOT.all_trans.(con_or_div).baseline>0);
% baseline_extant_trans=all_days_all_analysis.PLOT.all_trans.(con_or_div).baseline_extant_transitions; % for a shorter name


for i=1:num_of_days;
    syl_labels{i}=all_days_all_analysis.data{days_with_analysis(i)}.all_songs.all_trans.summary.syl_labels_in_order;
    all_days_all_analysis.PLOT.all_trans.(con_or_div).fractions{days_with_analysis(i)}=...
        all_days_all_analysis.data{days_with_analysis(i)}.all_songs.all_trans.summary.(con_or_div).matrix_of_fractions; 
    all_days_all_analysis.PLOT.all_trans.(con_or_div).fractions{days_with_analysis(i)}=100.*...
        all_days_all_analysis.PLOT.all_trans.(con_or_div).fractions{days_with_analysis(i)}; %convert from fraction to percent
    
%     get logicals telling whether this day had new or dropped transitions
%     relative to baseline)
%     all_days_all_analysis.PLOT.all_trans.(con_or_div).fractions{days_with_analysis(i)}(baseline_extant_trans)
    
end

   



% PLOT THE TRANSITION MATRICES
num_figures=ceil((num_of_days)/9); 
for i=1:num_figures;
    if i==num_figures;
        row_plots(i)=ceil((num_of_days-(i-1)*9)/3);
        col_plots(i)=ceil((num_of_days-(i-1)*9)/row_plots(i));
    else
        row_plots(i)=3;
        col_plots(i)=3;
    end
end

counter = 1; % count goes up 1 after each day plotted
while counter<=num_of_days;
    figure_num=ceil(counter/9);
    figure(figure_num); hold on;
    
    

    
    i=counter-(figure_num-1)*9; % index for specific figure
    j=counter; %index out of all days with data (same as i, but i resets with each figure)    
    syl_num=size(syl_labels{j},2);

        subplot(col_plots(figure_num),row_plots(figure_num),i);
        syl_mat = all_days_all_analysis.PLOT.all_trans.(con_or_div).fractions{days_with_analysis(j)};
%         max_value_local=max(max(all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre(:,:,days_with_analysis(j))));
%         min_value_local=min(min(all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre(:,:,days_with_analysis(j))));
        imagesc(syl_mat,[0 100]);
        colormap(flipud(gray))
        textStrings = num2str(syl_mat(:),'%0.0f');
        textStrings = strtrim(cellstr(textStrings));
        [x,y]=meshgrid(1:syl_num);
        hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
        midValue = mean(get(gca,'CLim'));
        textColors = repmat(syl_mat(:) > midValue,1,3);
        set(hStrings,{'Color'},num2cell(textColors,2));
        set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels{j}, 'YTick',1:syl_num,...
            'YTickLabel',syl_labels{j}, 'TickLength',[0 0]);
        ylabel('transition from:')
        xlabel('transition to:')
        title('\bf Baseline fractions (mean)');
        freezeColors

        set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels{j}, 'YTick',1:syl_num,...
            'YTickLabel',syl_labels{j}, 'TickLength',[0 0]);
        ylabel('transition from:')
        xlabel('transition to:')
        title(['day # ' num2str(days_with_analysis(j)) ', ' all_days_all_analysis.data{days_with_analysis(j)}.date]);
        freezeColors
        
        
        if counter == num_of_days;
            % make one colorbar
            ax=gca;
            pos=get(gca,'pos'); % to arrange colorbar relative to the last plot
            hc=colorbar('position',[pos(1)+1.2*pos(3) pos(2)+(0.2*pos(4)) 0.05*pos(3) 0.7*pos(4)]);
            % hc=colorbar('location','eastoutside')
            % set(hc,'xaxisloc','top');
            
            
            % give a general title
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,['\bf ' con_or_div ' transition matrices for each day'],'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top')
        end
        counter = counter+1;
end

% all_days_all_analysis.PLOT.all_trans.difference_matrices.common_syls=common_syls;
% all_days_all_analysis.PLOT.all_trans.difference_matrices.syl_indices=syl_indices;

% save
% cd PLOT/
% save('all_days_all_analysis',all_days_all_analysis

%% plot difference matrices for each day

% Figure out the syllables that were sang every day.
days_with_analysis=[];
for i=1:size(all_days_all_analysis.data,2);
    try
        if ~isempty(all_days_all_analysis.data{i}.all_songs.all_trans);
            days_with_analysis=[days_with_analysis i];
        end
    catch err
        continue
    end
end

clear syl_labels_concatenated
num_of_days=length(days_with_analysis);
for j=1:num_of_days;
    syl_labels_concatenated{j}=all_days_all_analysis.data{days_with_analysis(j)}.all_songs.all_trans.summary.syl_labels_in_order;
end

[common_syls syl_indices]=mintersect_array(syl_labels_concatenated); % find common syls and their indices in each day


% FIRST-part A, (con_or_div)
% get baseline
clear dummy
for j=1:length(baseline_days);
    dummy(:,:,j)=all_days_all_analysis.data{baseline_days(j)}.all_songs.all_trans.summary.(con_or_div).matrix_of_fractions(syl_indices{j},syl_indices{j});
end
all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline=mean(dummy,3);
all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline_extant_trans=...
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline>0;
all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline_extant_trans_sig=...
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline>0.05; % those that are over 5%

baseline_extant_trans=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline_extant_trans; % rename to make life easier
baseline_extant_trans_sig=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline_extant_trans_sig;

for i=1:num_of_days;
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre(:,:,days_with_analysis(i))=...
        all_days_all_analysis.data{days_with_analysis(i)}.all_songs.all_trans.summary.(con_or_div).matrix_of_fractions(syl_indices{i},syl_indices{i})...
        -all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline;
    
    % find transitions that are dropped or added on this day
    matrix=all_days_all_analysis.data{days_with_analysis(i)}.all_songs.all_trans.summary.(con_or_div).matrix_of_fractions(syl_indices{i},syl_indices{i});
    indicator=matrix>0;
    [add1 add2]=find(indicator>baseline_extant_trans);
    [drop1 drop2] = find(indicator<baseline_extant_trans);
    
    indicator_sig=matrix>=0.05; % find those new trans that are over 5%, or dropped trans that were >5%;
    [add1sig add2sig]=find(indicator_sig>baseline_extant_trans);
    [drop1sig drop2sig] = find(indicator<baseline_extant_trans_sig);
 
    
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).added_trns{days_with_analysis(i)}=[add1 add2];
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).dropped_trns{days_with_analysis(i)}=[drop1 drop2];
    
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).added_trns_sig{days_with_analysis(i)}=[add1sig add2sig];
    all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).dropped_trns_sig{days_with_analysis(i)}=[drop1sig drop2sig];

end

% convert from fraction to percent
all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre=... 
    100*all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre
all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline=...
    100*all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline;

% plot
clear row_plots
clear col_plots
num_figures2=ceil((num_of_days+1)/9); % +1 becuase have a baseline mean figure as well.
for i=1:num_figures2;
    if i==num_figures2;
        row_plots(i)=ceil((num_of_days+1-(i-1)*9)/3);
        col_plots(i)=ceil((num_of_days+1-(i-1)*9)/row_plots(i));
    else
        row_plots(i)=3;
        col_plots(i)=3;
    end
end

% 
% if num_of_days<6;
%     row_plots(1)=min(ceil(num_of_days+1)/2,4);
%     num_figures=1;
%     col_plots=ceil((num_of_days+1)/row_plots);
% elseif num_of_days>5 && num_of_days<9;
%     row_plots(1)=4;
%     num_figures=1;
%     col_plots=ceil((num_of_days+1)/row_plots);
% elseif num_of_days>8 && num_of_days<18;
%     num_figures=2;
%     row_plots(1)=4;
%     row_plots(2)=ceil((num_of_days-8)/3);
%     col_plots(1)=ceil((num_of_days+1)/row_plots(1));
%     col_plots(2)=ceil((num_of_days+1)/row_plots(2)); 
% end

counter = 1; % count goes up 1 after each day plotted
while counter<=num_of_days+1;
    figure_num=ceil(counter/9);
    figure_num2=figure_num+num_figures; % so I won't override the figures from 1st section.
    figure(figure_num2); hold on;
    
    syl_num=size(common_syls,2);
    syl_labels=common_syls;
    
    if counter==1; %only plot baseline for 1st figure
        %PLOT THE BASELINE FRACTION MATRIX
        subplot(col_plots(figure_num),row_plots(figure_num),1);
        syl_mat = all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).baseline;
        imagesc(syl_mat);
        colormap(flipud(gray))
        textStrings = num2str(syl_mat(:),'%0.0f');
        textStrings = strtrim(cellstr(textStrings));
        [x,y]=meshgrid(1:syl_num);
        hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
        midValue = mean(get(gca,'CLim'));
        textColors = repmat(syl_mat(:) > midValue,1,3);
        set(hStrings,{'Color'},num2cell(textColors,2));
        set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels, 'YTick',1:syl_num,...
            'YTickLabel',syl_labels, 'TickLength',[0 0]);
        ylabel('transition from:')
        xlabel('transition to:')
        title('\bf Baseline fractions (mean)');
        freezeColors
        counter=2;
    end
    
    % determine a common color scale across all days - normalize by max element
    % out of entire matrix
    max_value=max(max(max(all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre)));
    min_value=min(min(min(all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre)));
    absolute_max=max(min_value,max_value);
    
    % grayscale=[hot(100); flipud(cool(100))];
    % grayscale=cool(200)
    
    colorscale = [flipud([0.01:0.01:1; zeros(1,100); zeros(1,100)]'); [zeros(1,100); zeros(1,100); 0.01:0.01:1]'];
    
    %PLOT THE DIFFERENCE MATRICES
    i=counter-(figure_num-1)*9; % index for specific figure
    j=counter-1; %index out of all days with data
    subplot(col_plots(figure_num),row_plots(figure_num),i);
    syl_mat = all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre(:,:,days_with_analysis(j));
    max_value_local=max(max(all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre(:,:,days_with_analysis(j))));
    min_value_local=min(min(all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).all_days_incl_pre(:,:,days_with_analysis(j))));
    imagesc(syl_mat,[-absolute_max absolute_max]);
    %     colormap(grayscale(50+floor(min_value/2):50+ceil(max_value/2),:,:)); colorbar;
    colormap(colorscale);
    %     if i==num_of_days;
    %         h_bar=colorbar;
    %     end
    textStrings = num2str(syl_mat(:),'%0.0f');
    textStrings = strtrim(cellstr(textStrings));
    [x,y]=meshgrid(1:syl_num);
    hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
    %     midValue = mean(get(gca,'CLim'));
    textColors = repmat(syl_mat(:) >= 10,1,3)+repmat(syl_mat(:) <= -10,1,3); %only put in text of those that are past + or - 10.
    textColors_weak=0.35.*repmat(syl_mat(:) <10 & syl_mat(:) > -10,1,3);
    textColors=textColors+textColors_weak;
    
    set(hStrings,{'Color'},num2cell(textColors,2));
    %         set(hStrings,{'Color'},num2cell(textColors_weak,2));
    
    %     set(hStrings,{'Color'},{'w'})
    set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels, 'YTick',1:syl_num,...
        'YTickLabel',syl_labels, 'TickLength',[0 0]);
    ylabel('transition from:')
    xlabel('transition to:')
    title(['day # ' num2str(days_with_analysis(j)) ', ' all_days_all_analysis.data{days_with_analysis(j)}.date]);
    ax=gca;
    freezeColors
    
    % Annotate transitions that were added or dropped
    hold on; 
    xvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).added_trns{days_with_analysis(j)}(:,1);
    yvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).added_trns{days_with_analysis(j)}(:,2);
    scatter(yvals+0.3,xvals-0.3,'g');
    
    xvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).dropped_trns{days_with_analysis(j)}(:,1);
    yvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).dropped_trns{days_with_analysis(j)}(:,2);
    scatter(yvals+0.3,xvals+0.3,'r');
    
    %below, same but only sig drops and adds
    xvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).added_trns_sig{days_with_analysis(j)}(:,1);
    yvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).added_trns_sig{days_with_analysis(j)}(:,2);
    scatter(yvals+0.3,xvals-0.3,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    xvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).dropped_trns_sig{days_with_analysis(j)}(:,1);
    yvals=all_days_all_analysis.PLOT.all_trans.difference_matrices.(con_or_div).dropped_trns_sig{days_with_analysis(j)}(:,2);
    scatter(yvals+0.3,xvals+0.3,'MarkerEdgeColor','k','MarkerFaceColor','r');
        
    if counter == num_of_days+1;
        % make one colorbar
        ax=gca;
        pos=get(gca,'pos'); % to arrange colorbar relative to the last plot
        % set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);
        % pos=get(gca,'pos');
        % hc=colorbar('position',[pos(1) pos(2) pos(3) 0.5*pos4]);
        hc=colorbar('position',[pos(1)+1.2*pos(3) pos(2)+(0.2*pos(4)) 0.05*pos(3) 0.7*pos(4)]);
        % hc=colorbar('location','eastoutside')
        % set(hc,'xaxisloc','top');
        
        
        % give a general title
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,['\bf ' con_or_div ' transitions: difference from baseline (first ' num2str(max(baseline_days)) ' days)'],'HorizontalAlignment'...
            ,'center','VerticalAlignment', 'top')
    end

    counter = counter+1;
end

% all_days_all_analysis.PLOT.all_trans.difference_matrices.common_syls=common_syls;
% all_days_all_analysis.PLOT.all_trans.difference_matrices.syl_indices=syl_indices;

% save
% cd PLOT/
% save('all_days_all_analysis',all_days_all_analysis

%% Find new and dropped transitions and mark on plots

%% PLOT single divergent or convergent transition over days





end

