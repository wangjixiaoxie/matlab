%% LT 6/24/14 - Use with other lt_all_days_ContextSeq... code to analyse probe trials within a context seq learning experiment.
% Run this from another script. For example:
%         clear all;
%         load /bluejay3/lucas/birds/rd66gr93/all_days_transitions_ContextSeq_analysis/PLOT/31May2014_to_24Jun2014_ConAndDiv_ab_24Jun2014_1959/transitions_all_days_PLOT_31May2014_to_24Jun2014_ConAndDiv_ab_24Jun2014_1959
%
%         % INPUT PARAMETERS
%         % enter the probe epochs
%         probes_contextB = {'11Jun2014-1308_1556','16Jun2014-1235_1423','18Jun2014-1208_1332','19Jun2014-1219_1244','21Jun2014-1255_1519','24Jun2014-1239_1348'}; % these: WN off during CB
%         probes_contextA_early = {'23Jun2014-1138_1301'}; % these are WN played as probe in CA
%         % for plotting
%         edge_bin_size=6;
%         trn_to_plot='ab';
%
%         lt_all_days_ContextSeq_PROBE_CompilePlot


% Extract parameters
first_day=transitions_all_days.data{1}.parameters.date;
num_days=size(transitions_all_days.data,2);
last_day=transitions_all_days.data{num_days}.parameters.date;
[ProbeTimesCB]=lt_convert_EventTimes_to_RelTimes(first_day,probes_contextB); % extracts relevant information from the array containing probe times
[ProbeTimesCAe]=lt_convert_EventTimes_to_RelTimes(first_day,probes_contextA_early); % extracts relevant information from the array containing probe times
[ProbeTimes_WNoff_contextB_late]=lt_convert_EventTimes_to_RelTimes(first_day,probes_WNoff_contextB_late);
% NOTE: have not used ProbeTimes_WNoff_contextB_late to get data yet.  see
% directly below.

for i=1:num_days;
    % Extracting Context B probe trials
    if sum(ProbeTimesCB.JustDays_rel==i)==1; % if the current day (i) has probe data
        iii=find(ProbeTimesCB.JustDays_rel==i); % index within the probe array
        % first go from cell type probe time data to format in hours (e.g. 12:30pm = 12.50)
        probe_start_str=probes_contextB{iii}(1:14);
        probe_end_str=[probes_contextB{iii}(1:9) '-' probes_contextB{iii}(16:19)];
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_start_str,'ddmmmyyyy-HHMM'));
        ProbeStart=dummy.hours;
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_end_str,'ddmmmyyyy-HHMM'));
        ProbeEnd=dummy.hours;
        
        % find the epoch that is within those hours, pre, and post (all within same
        % context)
        ProbeIndsB_Dur{iii}=find(transitions_all_days.data{i}.contextB.time_hours>=ProbeStart & transitions_all_days.data{i}.contextB.time_hours<=ProbeEnd);
        ProbeIndsB_Post{iii}=find(transitions_all_days.data{i}.contextB.time_hours>ProbeEnd);
        % ADD exemptions for weird days
        if i==17; % today, only played one bout WN after probe, so no "post";
            ProbeIndsB_Post{iii}=[];
        end
        ProbeIndsB_Pre{iii}=find(transitions_all_days.data{i}.contextB.time_hours<ProbeStart);
        
        % extract all data, sorted for those contexts
        Fields=fieldnames(transitions_all_days.data{i}.contextB);
        for ii=1:length(Fields);
            try % use try since will get error if try to transfer over binned renditions
                ProbeData.ContextB.dur.(Fields{ii}){i} = transitions_all_days.data{i}.contextB.(Fields{ii})(ProbeIndsB_Dur{iii});
                ProbeData.ContextB.post.(Fields{ii}){i} = transitions_all_days.data{i}.contextB.(Fields{ii})(ProbeIndsB_Post{iii});
                ProbeData.ContextB.pre.(Fields{ii}){i} = transitions_all_days.data{i}.contextB.(Fields{ii})(ProbeIndsB_Pre{iii});
            catch err
                continue
            end
        end
    end
    
    
    % Extracting Context A probe trials
    if sum(ProbeTimesCAe.JustDays_rel==i)==1; % if the current day (i) has probe data
        iii=find(ProbeTimesCAe.JustDays_rel==i);
        % first go from cell type probe time data to format in hours (e.g. 12:30pm = 12.50)
        probe_start_str=probes_contextA_early{iii}(1:14);
        probe_end_str=[probes_contextA_early{iii}(1:9) '-' probes_contextA_early{iii}(16:19)];
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_start_str,'ddmmmyyyy-HHMM'));
        ProbeStart=dummy.hours;
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_end_str,'ddmmmyyyy-HHMM'));
        ProbeEnd=dummy.hours;
        
        % find the epoch that is within those hours, pre, and post (all within same
        % context)
        ProbeIndsAe_Dur{iii}=find(transitions_all_days.data{i}.contextA_early.time_hours>=ProbeStart & transitions_all_days.data{i}.contextA_early.time_hours<=ProbeEnd);
        ProbeIndsAe_Post{iii}=find(transitions_all_days.data{i}.contextA_early.time_hours>ProbeEnd);
        ProbeIndsAe_Pre{iii}=find(transitions_all_days.data{i}.contextA_early.time_hours<ProbeStart);
        
        % extract all data, sorted for those contexts
        Fields=fieldnames(transitions_all_days.data{i}.contextA_early);
        for ii=1:length(Fields);
            try % use try since will get error if try to transfer over binned renditions
                ProbeData.ContextA_early.dur.(Fields{ii}){i} = transitions_all_days.data{i}.contextA_early.(Fields{ii})(ProbeIndsAe_Dur{iii});
                ProbeData.ContextA_early.post.(Fields{ii}){i} = transitions_all_days.data{i}.contextA_early.(Fields{ii})(ProbeIndsAe_Post{iii});
                ProbeData.ContextA_early.pre.(Fields{ii}){i} = transitions_all_days.data{i}.contextA_early.(Fields{ii})(ProbeIndsAe_Pre{iii});
            catch err
                continue
            end
        end
    end
    end

% Assign stuff to output structure
if (1)
    try
        ProbeData.parameters.ProbeIndsB_Dur=ProbeIndsB_Dur;
        ProbeData.parameters.ProbeIndsB_Pre=ProbeIndsB_Pre;
        ProbeData.parameters.ProbeIndsB_Post=ProbeIndsB_Post;
    catch err
    end
    try
        ProbeData.parameters.ProbeIndsAe_Dur=ProbeIndsAe_Dur;
        ProbeData.parameters.ProbeIndsAe_Pre=ProbeIndsAe_Pre;
        ProbeData.parameters.ProbeIndsAe_Post=ProbeIndsAe_Post;
    catch err
    end
    
    ProbeData.parameters.Probedays_rel_CB=ProbeTimesCB.JustDays_rel;
    ProbeData.parameters.Probedays_rel_CAe=ProbeTimesCAe.JustDays_rel;
    ProbeData.parameters.probes_contextB=probes_contextB;
    ProbeData.parameters.probes_contextA_early=probes_contextA_early;
    ProbeData.parameters.probes_WNoff_contextB_late=probes_WNoff_contextB_late;

end


%% PLOT
% 1) baseline a to b
% 2) WN a to b
% 3) probe a to b
% 4) probe b to WN
% 5) false probe a to WN
% (not in that order necessarily)

contxt='ContextB';
% ProbePlotCols=lt_make_plot_colors(size(ProbeData.parameters.Probedays_rel_CB,1),1,[0 0 1]);

% PROBE DAYS (CB WN off)
figure; hold on; 
hfig1=subplot(1,5,1); xlim(Xlimits);ylim(Ylimits);hold on; title('CA-->CB')
hfig2=subplot(1,5,2); xlim(Xlimits);ylim(Ylimits);hold on; title('During probe')
hfig3=subplot(1,5,3); xlim(Xlimits);ylim(Ylimits);hold on; title('Probe-->WNon')
hfig3b=subplot(1,5,4); xlim(Xlimits);ylim(Ylimits);hold on;  title('During WN')
hfig3c=subplot(1,5,5); xlim(Xlimits);ylim(Ylimits);hold on;  title('WN-->CA')

Yall={[],[],[],[],[],[]};
Zall={[],[],[],[],[],[]};
tttmp_probe=[]; % added lt 4/13/15 to collect probe data for ttest
for i=1:size(ProbeData.parameters.Probedays_rel_CB,1)
    ii=ProbeData.parameters.Probedays_rel_CB(i);
    X={}; Y={}; Z={};
    % get numerator
    Y{ii}(1)=mean(transitions_all_days.data{ii}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
    if size(ProbeData.(contxt).dur.(trn_to_plot){ii},2)>=edge_bin_size;
        Y{ii}(2)=mean(ProbeData.(contxt).dur.(trn_to_plot){ii}(1:edge_bin_size)); % start of probe
        Y{ii}(3)=mean(ProbeData.(contxt).dur.(trn_to_plot){ii}(end-edge_bin_size+1:end)); % end of probe
        else
        Y{ii}(2)=mean(ProbeData.(contxt).dur.(trn_to_plot){ii}(1:end));
        Y{ii}(3)=mean(ProbeData.(contxt).dur.(trn_to_plot){ii}(1:end));
    end

    
    % to combine all renditions across days into one mean value
    Yall{1}=[Yall{1} transitions_all_days.data{ii}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end)];
        if size(ProbeData.(contxt).dur.(trn_to_plot){ii},2)>=edge_bin_size;
    Yall{2}=[Yall{2} ProbeData.(contxt).dur.(trn_to_plot){ii}(1:edge_bin_size)];
    Yall{3}=[Yall{3} ProbeData.(contxt).dur.(trn_to_plot){ii}(end-edge_bin_size+1:end)];
        else
    Yall{2}=[Yall{2} ProbeData.(contxt).dur.(trn_to_plot){ii}];
    Yall{3}=[Yall{3} ProbeData.(contxt).dur.(trn_to_plot){ii}];
        end

    % get demoninator
    Z{ii}(1)=mean(transitions_all_days.data{ii}.contextA_early.num_trans(end-edge_bin_size+1:end));
    if size(ProbeData.(contxt).dur.(trn_to_plot){ii},2)>=edge_bin_size;
        Z{ii}(2)=mean(ProbeData.(contxt).dur.num_trans{ii}(1:edge_bin_size));
        Z{ii}(3)=mean(ProbeData.(contxt).dur.num_trans{ii}(end-edge_bin_size+1:end));
    else
        Z{ii}(2)=mean(ProbeData.(contxt).dur.num_trans{ii});
        Z{ii}(3)=mean(ProbeData.(contxt).dur.num_trans{ii});
    end
    
    Zall{1}=[Zall{1} transitions_all_days.data{ii}.contextA_early.num_trans(end-edge_bin_size+1:end)];
    if size(ProbeData.(contxt).dur.(trn_to_plot){ii},2)>=edge_bin_size;
        Zall{2}=[Zall{2} ProbeData.(contxt).dur.num_trans{ii}(1:edge_bin_size)];
    Zall{3}=[Zall{3} ProbeData.(contxt).dur.num_trans{ii}(end-edge_bin_size+1:end)];
    else
             Zall{2}=[Zall{2} ProbeData.(contxt).dur.num_trans{ii}];
    Zall{3}=[Zall{3} ProbeData.(contxt).dur.num_trans{ii}];
    end
        X{ii}=Y{ii}./Z{ii};

    
        try % try becuase sometimes does not have post (i.e. kep probe on entire context)
    Y{ii}(4)=mean(ProbeData.(contxt).post.(trn_to_plot){ii}(1:edge_bin_size)); % after WN start
    Y{ii}(5)=mean(ProbeData.(contxt).post.(trn_to_plot){ii}(end-edge_bin_size+1:end)); % end of WN (i.e. CB)
    
    Yall{4}=[Yall{4} ProbeData.(contxt).post.(trn_to_plot){ii}(1:edge_bin_size)];
    Yall{5}=[Yall{5} ProbeData.(contxt).post.(trn_to_plot){ii}(end-edge_bin_size+1:end)];
    
    Z{ii}(4)=mean(ProbeData.(contxt).post.num_trans{ii}(1:edge_bin_size));
    Z{ii}(5)=mean(ProbeData.(contxt).post.num_trans{ii}(end-edge_bin_size+1:end)); % end of WN (i.e. CB)
    
    
    Zall{4}=[Zall{4} ProbeData.(contxt).post.num_trans{ii}(1:edge_bin_size)];
    Zall{5}=[Zall{5} ProbeData.(contxt).post.num_trans{ii}(end-edge_bin_size+1:end)];
    
    
    
catch err
        end
try
    
        Y{ii}(6)=mean(transitions_all_days.data{ii}.contextA_late.(trn_to_plot)(1:edge_bin_size));  % beginning of CA
    Yall{6}=[Yall{6} transitions_all_days.data{ii}.contextA_late.(trn_to_plot)(1:edge_bin_size)];

         Z{ii}(6)=mean(transitions_all_days.data{ii}.contextA_late.num_trans(1:edge_bin_size));  % beginning of CA
       Zall{6}=[Zall{6} transitions_all_days.data{ii}.contextA_late.num_trans(1:edge_bin_size)];
catch err
end

       X{ii}=Y{ii}./Z{ii};

    
    % plot
    date=transitions_all_days.data{ii}.parameters.date(1:5);
    subplot(hfig1); plot([2 3],X{ii}(1:2),'-ob'); text(3.05,X{ii}(2),date);
    tttmp_probe=[tttmp_probe X{ii}(2)-X{ii}(1)];
    subplot(hfig2); plot([2 3],X{ii}(2:3),'-ob'); text(3.05,X{ii}(3),date);
    try 
        subplot(hfig3); plot([2 3],X{ii}(3:4),'-or'); text(3.05,X{ii}(4),date);
        subplot(hfig3b); plot([2 3],X{ii}(4:5),'-or'); text(3.05,X{ii}(5),date);
        subplot(hfig3c); plot([2 3],X{ii}(5:6),'-og'); text(3.05,X{ii}(6),date);
    catch err
    end
    
       try 
        subplot(hfig3c); plot([2 3],X{ii}(5:6),'-og'); text(3.05,X{ii}(6),date);
    catch err
    end
 
end
% compile across all days (both across renditions, or keeping day info intact)
Xall=[sum(Yall{1})/sum(Zall{1}) sum(Yall{2})/sum(Zall{2}) sum(Yall{3})/sum(Zall{3}) sum(Yall{4})/sum(Zall{4})...
    sum(Yall{5})/sum(Zall{5}) sum(Yall{6})/sum(Zall{6})];
subplot(hfig1); plot([2 3],Xall(1:2),'-ob','MarkerFaceColor','b');
subplot(hfig2); plot([2 3],Xall(2:3),'-ob','MarkerFaceColor','b');
subplot(hfig3); plot([2 3],Xall(3:4),'-or','MarkerFaceColor','r');
subplot(hfig3b); plot([2 3],Xall(4:5),'-or','MarkerFaceColor','r');
subplot(hfig3c); plot([2 3],Xall(5:6),'-og','MarkerFaceColor','g');


% ADDING global title
[ax, h]=subtitle(['WN off PROBE (CB); Transition probability of ' trn_to_plot ' (' transitions_all_days.parameters.phrase ') at edges, bin size: '  num2str(edge_bin_size) ' songs']);
set(h,'Position',title_pos);



% FALSE PROBE DAYS (i.e. WN in context A)
% NOTE: NEED TO FIX. look into whether I should have both beginning and end of CA (postprobe)
figure; 
    hfig4=subplot(1,7,1); xlim(Xlimits); ylim(Ylimits);hold on; title('WNoff-->WNon') % false probes context A
    hfig5=subplot(1,7,2); xlim(Xlimits); ylim(Ylimits);hold on; title('During probe') % false probes context A
    hfig6=subplot(1,7,3); xlim(Xlimits); ylim(Ylimits);hold on; title('Probe-->WNoff (in CA)') % false probes context A
    hfig6b=subplot(1,7,4); xlim(Xlimits); ylim(Ylimits);hold on; title('Probe-->CB immedate (if happened)') %
    hfig6c=subplot(1,7,5); xlim(Xlimits); ylim(Ylimits);hold on; title('CA(non-probe)-->CB (if happened)') %
   
    hfig6d=subplot(1,7,6); xlim(Xlimits); ylim(Ylimits);hold on; title('During CB') %
    hfig6e=subplot(1,7,7); xlim(Xlimits); ylim(Ylimits);hold on; title('CB-->CA') %
    
    Yall={[],[],[],[],[],[],[]};
    Zall={[],[],[],[],[],[],[]};
    
    for i=1:size(ProbeData.parameters.Probedays_rel_CAe,1)
        ii=ProbeData.parameters.Probedays_rel_CAe(i);
        X=[]; Y=[]; Z=[];
        % get numerator
        Y(1)=mean(ProbeData.ContextA_early.pre.(trn_to_plot){ii}(end-edge_bin_size+1:end)); % pre
        Y(2)=mean(ProbeData.ContextA_early.dur.(trn_to_plot){ii}(1:edge_bin_size)); % dur(start)
        Y(3)=mean(ProbeData.ContextA_early.dur.(trn_to_plot){ii}(end-edge_bin_size+1:end)); % dur(end)
        
        Yall{1}=[Yall{1} ProbeData.ContextA_early.pre.(trn_to_plot){ii}(end-edge_bin_size+1:end)];
        Yall{2}=[Yall{2} ProbeData.ContextA_early.dur.(trn_to_plot){ii}(1:edge_bin_size)];
        Yall{3}=[Yall{3} ProbeData.ContextA_early.dur.(trn_to_plot){ii}(end-edge_bin_size+1:end)];
        
        % get demoninator
        Z(1)=mean(ProbeData.ContextA_early.pre.num_trans{ii}(end-edge_bin_size+1:end)); % pre
        Z(2)=mean(ProbeData.ContextA_early.dur.num_trans{ii}(1:edge_bin_size)); % dur(start)
        Z(3)=mean(ProbeData.ContextA_early.dur.num_trans{ii}(end-edge_bin_size+1:end)); % dur(end)
        
        Zall{1}=[Zall{1} ProbeData.ContextA_early.pre.num_trans{ii}(end-edge_bin_size+1:end)];
        Zall{2}=[Zall{2} ProbeData.ContextA_early.dur.num_trans{ii}(1:edge_bin_size)];
        Zall{3}=[Zall{3} ProbeData.ContextA_early.dur.num_trans{ii}(end-edge_bin_size+1:end)];
        
        X=Y./Z;
        
        % plot
        date=transitions_all_days.data{ii}.parameters.date(1:5);
        subplot(hfig4); plot([2 3],X(1:2),'-or'); text(3.05,X(2),date);
        subplot(hfig5); plot([2 3],X(2:3),'-or'); text(3.05,X(3),date);
        
        try % try plotting end of probe to post (sometimes lacks - goes straight to context b)
            Y(4)=mean(ProbeData.ContextA_early.post.(trn_to_plot){ii}(1:edge_bin_size)); % post (still in CA.  if move to CB then blank)
            Yall{4}=[Yall{4} ProbeData.ContextA_early.post.(trn_to_plot){ii}(1:edge_bin_size)];
            Z(4)=mean(ProbeData.ContextA_early.post.num_trans{ii}(1:edge_bin_size)); % post (still in CA.  if move to CB then blank)
            Zall{4}=[Zall{4} ProbeData.ContextA_early.post.num_trans{ii}(1:edge_bin_size)];
            X=Y./Z;
            
        catch err
            
        end
        
        try % try plotting end of probe straight to CB
            Y(5)=mean(transitions_all_days.data{ii}.contextB.(trn_to_plot)(1:edge_bin_size)); % start CB
            Yall{5}=[Yall{5} transitions_all_days.data{ii}.contextB.(trn_to_plot)(1:edge_bin_size)];
            
            Z(5)=mean(transitions_all_days.data{ii}.contextB.num_trans(1:edge_bin_size)); % post (still in CA.  if move to CB then blank)
            Zall{5}=[Zall{5} transitions_all_days.data{ii}.contextB.num_trans(1:edge_bin_size)];
            X=Y./Z;
            
            if X(4)>0; % i.e. if there is any song between probe and CB
                subplot(hfig6); plot([2 3],X(3:4),'-og'); text(3.05,X(4),date); % end of probe to CA
                subplot(hfig6c); plot([2 3],X(4:5),'-or'); text(3.05,X(5),date); % CA to CB
            else
                subplot(hfig6b); plot([2 3],X([3 5]),'-or'); text(3.05,X(5),date);
            end
        catch err
        end
        
        try % try plotting during CB and Cb to CA
            Y(6)=mean(transitions_all_days.data{ii}.contextB.(trn_to_plot)(end-edge_bin_size+1:end)); % end of CB
            Y(7)=mean(transitions_all_days.data{ii}.contextA_late.(trn_to_plot)(1:edge_bin_size)); % start CAlate
            
            
            Yall{6}=[Yall{6} transitions_all_days.data{ii}.contextB.(trn_to_plot)(end-edge_bin_size+1:end)];
            Yall{7}=[Yall{7} transitions_all_days.data{ii}.contextA_late.(trn_to_plot)(1:edge_bin_size)];
            
            
            Z(6)=mean(transitions_all_days.data{ii}.contextB.num_trans(end-edge_bin_size+1:end));
            Z(7)=mean(transitions_all_days.data{ii}.contextA_late.num_trans(1:edge_bin_size));
            
            Zall{6}=[Zall{6} transitions_all_days.data{ii}.contextB.num_trans(end-edge_bin_size+1:end)];
            Zall{7}=[Zall{7} transitions_all_days.data{ii}.contextA_late.num_trans(1:edge_bin_size)];
            
            X=Y./Z;
            
            subplot(hfig6d); plot([2 3],X(5:6),'-or'); text(3.05,X(6),date);
            subplot(hfig6e); plot([2 3],X(6:7),'-og'); text(3.05,X(7),date);
            
        catch err
        end
        
    end
    
    
    Xall=[sum(Yall{1})/sum(Zall{1}) sum(Yall{2})/sum(Zall{2}) sum(Yall{3})/sum(Zall{3}) sum(Yall{4})/sum(Zall{4})...
        sum(Yall{5})/sum(Zall{5}) sum(Yall{6})/sum(Zall{6}) sum(Yall{7})/sum(Zall{7})];
    subplot(hfig4); plot([2 3],Xall(1:2),'-or','MarkerFaceColor','r');
    subplot(hfig5); plot([2 3],Xall(2:3),'-or','MarkerFaceColor','r');
    subplot(hfig6); plot([2 3],Xall(3:4),'-og','MarkerFaceColor','g');
    
    subplot(hfig6c); plot([2 3],Xall(4:5),'-or','MarkerFaceColor','r');
%     subplot(hfig6b); plot([2 3],Xall([3 5]),'-or','MarkerFaceColor','r'); modify so that the left dot does not contain all days.
    
    subplot(hfig6d); plot([2 3],Xall(5:6),'-or','MarkerFaceColor','r');
    subplot(hfig6e); plot([2 3],Xall(6:7),'-og','MarkerFaceColor','r');
   
[ax, h]=subtitle(['WN on PROBE (CA); Transition probability of ' trn_to_plot ' (' transitions_all_days.parameters.phrase ') at edges, bin size: '  num2str(edge_bin_size) ' songs']);
set(h,'Position',title_pos);

    
    
    
% WN and baseline days
figure; 
% baseline days
hfig7=subplot(1,4,1); xlim(Xlimits); ylim(Ylimits);hold on; title('Baseline: CA-->CB');
Yall={[],[]};
Zall={[],[]};
for i=transitions_all_days.PLOT.parameters.BaselineDays;
    X=[]; Y=[]; Z=[];
    try
        
    if strcmp(contxt,'ContextB')==1; % then relevant baseline is CA to CB transition
        Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
        Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
        
        Yall{1}=[Yall{1} transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end)];
        Yall{2}=[Yall{2} transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size)];
        
        
        % get demoninator
        Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
        Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));
        Zall{1}=[Zall{1} transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end)];
        Zall{2}=[Zall{2} transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size)];
        
        X=Y./Z;
        
        date=transitions_all_days.data{i}.parameters.date(1:5);
        subplot(hfig7); plot([2 3],X,'-ok'); text(3.05,X(2),date);
        
    end
    catch err % not enough songs
    end
end
Xall=[sum(Yall{1})/sum(Zall{1}) sum(Yall{2})/sum(Zall{2})];
subplot(hfig7); plot([2 3],Xall,'-ok','MarkerFaceColor','k');



% WN days
hfig8=subplot(1,4,2); xlim(Xlimits); ylim(Ylimits);hold on; title('Training: CA-->CB');
hfig9=subplot(1,4,3); xlim(Xlimits); ylim(Ylimits);hold on;title('Training: During CB');
hfig10=subplot(1,4,4); xlim(Xlimits); ylim(Ylimits);hold on;title('Training: CB-->CA');
Yall={[],[],[],[]};
Zall={[],[],[],[]};

for i=transitions_all_days.PLOT.parameters.WNDays;
    if i~=11 && i~=13 ; % avoid these days becuase the 11) lights came on late 13) lights did not go out prev night.
        X=[]; Y=[]; Z=[];
        try  % not enough songs
            if isfield(transitions_all_days.data{i}, 'contextB')
                Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
                Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
                Y(3)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(end-edge_bin_size+1:end));
                Yall{1}=[Yall{1} transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end)];
                Yall{2}=[Yall{2} transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size)];
                Yall{3}=[Yall{3} transitions_all_days.data{i}.contextB.(trn_to_plot)(end-edge_bin_size+1:end)];
                try
                    Y(4)=mean(transitions_all_days.data{i}.contextA_late.(trn_to_plot)(1:edge_bin_size));
                    Yall{4}=[Yall{4} transitions_all_days.data{i}.contextA_late.(trn_to_plot)(1:edge_bin_size)];
                catch err % probably have not labeled data
                    Y(4)=nan;
                end
                
                % get demoninator
                Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
                Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));
                Z(3)=mean(transitions_all_days.data{i}.contextB.num_trans(end-edge_bin_size+1:end));
                
                Zall{1}=[Zall{1} transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end)];
                Zall{2}=[Zall{2} transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size)];
                Zall{3}=[Zall{3} transitions_all_days.data{i}.contextB.num_trans(end-edge_bin_size+1:end)];
                
                try
                    Z(4)=mean(transitions_all_days.data{i}.contextA_late.num_trans(1:edge_bin_size));
                    Zall{4}=[Zall{4} transitions_all_days.data{i}.contextA_late.num_trans(1:edge_bin_size)];
                catch err
                    Z(4)=nan;
                end
                
                X=Y./Z;
                
                % plot
                date=transitions_all_days.data{i}.parameters.date(1:5);
                subplot(hfig8); plot([2 3],X(1:2),'-or','MarkerFaceColor','r'); text(3.05,X(2),date);
                subplot(hfig9); plot([2 3],X(2:3),'-or','MarkerFaceColor','r'); text(3.05,X(3),date);
                subplot(hfig10); plot([2 3],X(3:4),'-og','MarkerFaceColor','g'); text(3.05,X(4),date);
            end
        catch err
        end
    end
end

Xall=[sum(Yall{1})/sum(Zall{1}) sum(Yall{2})/sum(Zall{2}) sum(Yall{3})/sum(Zall{3}) sum(Yall{4})/sum(Zall{4})];
subplot(hfig8); plot([2 3],Xall(1:2),'-ok','MarkerFaceColor','k');
subplot(hfig9); plot([2 3],Xall(2:3),'-ok','MarkerFaceColor','k');
subplot(hfig10); plot([2 3],Xall(3:4),'-ok','MarkerFaceColor','k');

% ADDING global title
[ax, h]=subtitle(['Non-probe days: Transition probability of ' trn_to_plot ' (' transitions_all_days.parameters.phrase ') at edges, bin size: '  num2str(edge_bin_size) ' songs']);
set(h,'Position',title_pos);


%% PLOTTING OVER DAYS THE DIFFERENCES
% make one 2 x (days) array that has pre and post context movement for all
% days (i.e. baseline, normal, and probe days)
X_all_days=[];
% Baseline
for i=transitions_all_days.PLOT.parameters.BaselineDays;
    X=[]; Y=[]; Z=[];
    try
        
    if strcmp(contxt,'ContextB')==1; % then relevant baseline is CA to CB transition
        
        Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
        Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
        
        % get demoninator
        Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
        Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));

        X=Y./Z;
X_all_days(:,i)=X';
        
    end
    catch err % not enough songs
    end
end

% WN days (non-probe)
for i=transitions_all_days.PLOT.parameters.WNDays;
%     if i~=11 && i~=13 ; % avoid these days becuase the 11) lights came on late 13) lights did not go out prev night.
        X=[]; Y=[]; Z=[];
        try  % not enough songs
            if isfield(transitions_all_days.data{i}, 'contextB')
                Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
                Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
                
                
                % get demoninator
                Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
                Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));
                
                                
                X=Y./Z;
            X_all_days(:,i)=X';
                
            end

        catch err
        end
%     end
end

% WN days (within probes, but still usable, since probe is not in A->B
% transition)
for ii=1:size(ProbeTimes_WNoff_contextB_late.JustDays_rel,1);
    i=ProbeTimes_WNoff_contextB_late.JustDays_rel(ii);
    X=[]; Y=[]; Z=[];
    try  % not enough songs
        if isfield(transitions_all_days.data{i}, 'contextB')
            Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
            Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
            
            % get demoninator
            Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
            Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));
            
            X=Y./Z;
            X_all_days(:,i)=X';
        end
    catch err
    end
end

% % WN days (Probe A days that went back to WN before Context B    
%     for i=1:size(ProbeData.parameters.Probedays_rel_CAe,1);
%         ii=ProbeData.parameters.Probedays_rel_CAe(i);
%         X=[]; Y=[]; Z=[];
%         % get numerator
%         try 
%             % try plotting end of probe to post (sometimes lacks - goes straight to context b)
%             Y(4)=mean(ProbeData.ContextA_early.post.(trn_to_plot){ii}(1:edge_bin_size)); % post (still in CA.  if move to CB then blank)
%             Z(4)=mean(ProbeData.ContextA_early.post.num_trans{ii}(1:edge_bin_size)); % post (still in CA.  if move to CB then blank)
%         % try plotting end of probe straight to CB
%             Y(5)=mean(transitions_all_days.data{ii}.contextB.(trn_to_plot)(1:edge_bin_size)); % start CB
%             Z(5)=mean(transitions_all_days.data{ii}.contextB.num_trans(1:edge_bin_size)); % post (still in CA.  if move to CB then blank)
%             X(1)=Y(4)/Z(4);
%             X(2)=Y(5)/Z(5);
%             
%             X_all_days(:,ii)=X';
%         catch err
%         end
%     end
    
    
% Probe days (WN off early B)
for i=1:size(ProbeData.parameters.Probedays_rel_CB,1)
    ii=ProbeData.parameters.Probedays_rel_CB(i);
    X=[]; Y=[]; Z=[];
    try
        % get numerator
        Y(1)=mean(transitions_all_days.data{ii}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
        if size(ProbeData.(contxt).dur.(trn_to_plot){ii},2)>=edge_bin_size;
            Y(2)=mean(ProbeData.(contxt).dur.(trn_to_plot){ii}(1:edge_bin_size)); % start of probe
        else
            Y(2)=mean(ProbeData.(contxt).dur.(trn_to_plot){ii}(1:end));
        end
        
        % get demoninator
        Z(1)=mean(transitions_all_days.data{ii}.contextA_early.num_trans(end-edge_bin_size+1:end));
        if size(ProbeData.(contxt).dur.(trn_to_plot){ii},2)>=edge_bin_size;
            Z(2)=mean(ProbeData.(contxt).dur.num_trans{ii}(1:edge_bin_size));
        else
            Z(2)=mean(ProbeData.(contxt).dur.num_trans{ii});
        end
        
        X=Y./Z;
        X_all_days(:,ii)=X';
    catch err
    end
end


% SECOND TRIAL - baseline
for i=transitions_all_days.PLOT.parameters.BaslineDays_trial2;
        X=[]; Y=[]; Z=[];
        try  % not enough songs
            if isfield(transitions_all_days.data{i}, 'contextB')
                Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
                Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
                
                
                % get demoninator
                Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
                Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));
                
                                
                X=Y./Z;
            X_all_days(:,i)=X';
                
            end

        catch err
        end
end

% SECOND TRIAL - WN
for i=transitions_all_days.PLOT.parameters.WNDays_trial2;
        X=[]; Y=[]; Z=[];
        try  % not enough songs
            if isfield(transitions_all_days.data{i}, 'contextB')
                Y(1)=mean(transitions_all_days.data{i}.contextA_early.(trn_to_plot)(end-edge_bin_size+1:end));
                Y(2)=mean(transitions_all_days.data{i}.contextB.(trn_to_plot)(1:edge_bin_size));
                
                
                % get demoninator
                Z(1)=mean(transitions_all_days.data{i}.contextA_early.num_trans(end-edge_bin_size+1:end));
                Z(2)=mean(transitions_all_days.data{i}.contextB.num_trans(1:edge_bin_size));
                
                                
                X=Y./Z;
            X_all_days(:,i)=X';
                
            end

        catch err
        end
end


% PLOT ALL
figure; hold on;
BL_days=transitions_all_days.PLOT.parameters.BaselineDays;
Training_days=transitions_all_days.PLOT.parameters.WNDays;
Probe_days_BWNoffearly=ProbeData.parameters.Probedays_rel_CB;




for i=1:length(X_all_days);
    try
        if find(Probe_days_BWNoffearly==i)>0;
        plot([i-0.3 i+0.3],X_all_days(:,i),'>-','Color','b','MarkerFaceColor','b','MarkerSize',6);
        elseif find(BL_days==i)>0;
        plot([i-0.3 i+0.3],X_all_days(:,i),'>-','Color','k','MarkerFaceColor','k','MarkerSize',6);        
        elseif find(ProbeData.parameters.Probedays_rel_CAe==i)>0;
        plot([i-0.3 i+0.3],X_all_days(:,i),'>-','Color','k','MarkerFaceColor','y','MarkerSize',6);
        elseif find(transitions_all_days.PLOT.parameters.BaslineDays_trial2==i)>0;
        plot([i-0.3 i+0.3],X_all_days(:,i),'>-','Color','k','MarkerFaceColor','k','MarkerSize',6);
        elseif find(transitions_all_days.PLOT.parameters.WNDays_trial2==i)>0;
        plot([i-0.3 i+0.3],X_all_days(:,i),'>-','Color','r','MarkerFaceColor','r','MarkerSize',6);
        else
        plot([i-0.3 i+0.3],X_all_days(:,i),'>-','Color','r','MarkerFaceColor','r','MarkerSize',6);
%         plot([i+0.3],X_all_days(2,i),'<-','Color','k','MarkerFaceColor','r');
       
        end      
    catch err
    end
end

% baseline line
line([BL_days(end)+0.5 BL_days(end)+0.5], ylim,'Color','k');
line([transitions_all_days.PLOT.parameters.BaslineDays_trial2(1)-0.5 transitions_all_days.PLOT.parameters.BaslineDays_trial2(1)-0.5], ylim,'Color','k');
line([transitions_all_days.PLOT.parameters.BaslineDays_trial2(end)+0.5 transitions_all_days.PLOT.parameters.BaslineDays_trial2(end)+0.5], ylim,'Color','k');

xlabel('days')
ylabel('transition probability');



%% SAVE
if save_yes_no==1;
cd(savedir);
try 
    cd ProbeData;
catch err
    mkdir ProbeData;
    cd ProbeData;
end

ProbeData.parameters.edge_bin_size=edge_bin_size;
ProbeData.parameters.trn_to_plot=trn_to_plot;
ProbeData.parameters.savedir=savedir;
ProbeData.parameters.filename=filename;

transitions_all_days.ProbeData=ProbeData;
save_label2=['_' save_label];

save(['transitions_all_days_ProbeData' save_label2],'transitions_all_days');
savemultfigs
end