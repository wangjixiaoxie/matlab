%% For context B

% CHANGE THIS NUMBER
ProbeTimesCAe.JustDays_rel % this gives you indices that you need to go thru manially
i=20 % index

% PROBE TIMES
        iii=find(ProbeTimesCB.JustDays_rel==i); % index within the probe array
        % first go from cell type probe time data to format in hours (e.g. 12:30pm = 12.50)
        probe_start_str=probes_contextB{iii}(1:14);
        probe_end_str=[probes_contextB{iii}(1:9) '-' probes_contextB{iii}(16:19)];
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_start_str,'ddmmmyyyy-HHMM'));
        ProbeStart=dummy.hours;
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_end_str,'ddmmmyyyy-HHMM'));
        ProbeEnd=dummy.hours;


        % PROBE song indices
        find(transitions_all_days.data{i}.contextB.time_hours>=ProbeStart & transitions_all_days.data{i}.contextB.time_hours<=ProbeEnd)
        
% date
transitions_all_days.data{i}.parameters.date

% things to visualize

transitions_all_days.data{i}.contextB.time_hours

ProbeStart

ProbeEnd

ProbeIndsB_Dur{iii}
ProbeIndsB_Post{iii}


% time conversions
.8053 * 60


%% for Context A

% CHANGE THIS NUMBER
ProbeTimesCAe.JustDays_rel % this gives you indices that you need to go thru manially
i=26 % index

% PROBE TIMES
        iii=find(ProbeTimesCAe.JustDays_rel==i);
        % first go from cell type probe time data to format in hours (e.g. 12:30pm = 12.50)
        probe_start_str=probes_contextA_early{iii}(1:14);
        probe_end_str=[probes_contextA_early{iii}(1:9) '-' probes_contextA_early{iii}(16:19)];
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_start_str,'ddmmmyyyy-HHMM'));
        ProbeStart=dummy.hours;
        [~, dummy]=lt_convert_datenum_to_hour(datenum(probe_end_str,'ddmmmyyyy-HHMM'));
        ProbeEnd=dummy.hours;


        % PROBE song indices
        find(transitions_all_days.data{i}.contextA_early.time_hours>=ProbeStart & transitions_all_days.data{i}.contextA_early.time_hours<=ProbeEnd);
        
% date
transitions_all_days.data{i}.parameters.date

% things to visualize

transitions_all_days.data{i}.contextA_early.time_hours

ProbeStart

ProbeEnd

ProbeIndsAe_Pre{iii}
ProbeIndsAe_Dur{iii}
ProbeIndsAe_Post{iii}


% time conversions
.8053 * 60
