function compiled_seqdep_pitch=lt_compile_seq_dep_pitch_data_SEQFILTER(compiled_seqdep_pitch,SeqPreList,SylTargList,SeqPostList)
%% LT 11/13/14 - Use on structure made with function lt_compile_seq_dep_pitch_data, to filter syls(e.g. b) to sequence specific syl (e.g. a[b]);

% INPUTS
% compiled_seqdep_pitch - structure that is output from
% lt_compile_seq_dep_pitch_data. This contains data for the single syls
% (e.g. b) of pitch, contour, sequence, etc.
% SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should
% combine. in this case the two sequences are ac[c]b and ab[b].
% SylTargList={'c','b'}; % these must already have raw data compiled above
% SeqPostList={'b',''};
%
% OUTPUTS:
% compiled_seqdep_pitch - overrides input structure by adding fields
% related tot the new sequences

curr_dir=pwd;

%% FILTER syllables based on sequence context (e.g. only b from ab);
for j=1:length(SylTargList);
    SeqPre=SeqPreList{j};
    SylTarg=SylTargList{j};
    SeqPost=SeqPostList{j};
    
    % Automatic params
    TargStr=[SeqPre SylTarg SeqPost]; % entire string
    TartPos=length(SeqPre)+1; % target syl position in string
    
    % field names for structure - capitalize the targ syl only
    SeqFieldName=[SeqPre upper(SylTarg) SeqPost];
    
    % First make sure this targ syl is present this day.
    
    
    % DATA WITH OUTLIERS REMOVED
    NumRends=size(compiled_seqdep_pitch.data.(SylTarg),1);
    IndsToKeep=[];
    for ii=1:NumRends; % all rends of the targ syl
        sylpos=compiled_seqdep_pitch.data.(SylTarg){ii,8}; % ind in song of candidate
        IndsToCheck=(sylpos-length(SeqPre)):(sylpos+length(SeqPost)); % what within song inds to look at
        
        if IndsToCheck(1)>0 && IndsToCheck(end)<=length(compiled_seqdep_pitch.data.(SylTarg){ii,7}); % sequence is expected to be within the 1st and last indices of the actual data.
            SeqToCheck=compiled_seqdep_pitch.data.(SylTarg){ii,7}(IndsToCheck); % is this the desired string?
            
            if strcmp(TargStr,SeqToCheck)==1; % if is desired string, then keep this datapt
                IndsToKeep=[IndsToKeep ii];
            end
        end
    end
    
    if ~isempty(IndsToKeep);
        compiled_seqdep_pitch.data.(SeqFieldName)=...
            compiled_seqdep_pitch.data.(SylTarg)(IndsToKeep,:);
    end        
    
    % DECIDED NOT to filter data that has not had outliers removed. - to
    % save file space
    
    % OUTLIERS REMOVED
    %     NumRends=size(compiled_seqdep_pitch.data_WithOutlier.(SylTarg),1);
    %     IndsToKeep=[];
    %     for ii=1:NumRends;
    %         sylpos=compiled_seqdep_pitch.data_WithOutlier.(SylTarg){ii,8}; % ind in song of candidate
    %         IndsToCheck=(sylpos-length(SeqPre)):(sylpos+length(SeqPost)); % what within song inds to look at
    %
    %         if IndsToCheck(1)>0 && IndsToCheck(end)<=length(compiled_seqdep_pitch.data.(SylTarg){ii,7});
    %             SeqToCheck=compiled_seqdep_pitch.data_WithOutlier.(SylTarg){ii,7}(IndsToCheck); % is this the desired string?
    %
    %             if strcmp(TargStr,SeqToCheck)==1; % if is desired string, then keep this datapt
    %                 IndsToKeep=[IndsToKeep ii];
    %             end
    %         end
    %     end
    %
    %     compiled_seqdep_pitch.data_WithOutlier.(SeqFieldName)=...
    %         compiled_seqdep_pitch.data_WithOutlier.(SylTarg)(IndsToKeep,:);
end




%% SUMMARY STATISTICS - e.g. mean, SD, etc.


syl_fields=fields(compiled_seqdep_pitch.data); % what seqs are there?
for i=1:length(syl_fields);
    
    % OUTLIERS REMOVED data.
    X=cell2mat(compiled_seqdep_pitch.data.(syl_fields{i})(:,1));
    compiled_seqdep_pitch.summary_stats.(syl_fields{i}).meanFF=mean(X);
    compiled_seqdep_pitch.summary_stats.(syl_fields{i}).medianFF=median(X);
    compiled_seqdep_pitch.summary_stats.(syl_fields{i}).sdFF=std(X);
    compiled_seqdep_pitch.summary_stats.(syl_fields{i}).n=length(X);
    compiled_seqdep_pitch.summary_stats.(syl_fields{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(compiled_seqdep_pitch.data.(syl_fields{i})(:,2)); % rows: rend, cols: timebin
    compiled_seqdep_pitch.summary_stats.(syl_fields{i}).pitchcountour_mean=mean(X,1);
    
    %     WILL IGNORE GETTING MEAN SPECTROGRAM - large file size.
    % get mean spectrogram
    %     X=[];
    %     for ii=1:size(compiled_seqdep_pitch.data.(syl_fields{i})(:,4),1); % num rends
    %         X(:,:,ii)=compiled_seqdep_pitch.data.(syl_fields{i}){ii,4}; % extract individual specs
    %     end
    %     compiled_seqdep_pitch.summary_stats.(syl_fields{i}).spec_mean=mean(X,3);
    
    
    % IGNORE DATA THAT HAS OUTLIERS.
    %     % NO OUTLIERS
    %     X=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syl_fields{i})(:,1));
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).meanFF=mean(X);
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).medianFF=median(X);
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).sdFF=std(X);
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).n=length(X);
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).semFF=std(X)/sqrt(length(X)-1);
    %
    %     % get mean pitch contour
    %     X=cell2mat(compiled_seqdep_pitch.data_WithOutlier.(syl_fields{i})(:,2)); % rows: rend, cols: timebin
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).pitchcountour_mean=mean(X,1);
    %
    %     % get mean spectrogram
    %     X=[];
    %     for ii=1:size(compiled_seqdep_pitch.data_WithOutlier.(syl_fields{i})(:,4),1); % num rends
    %         X(:,:,ii)=compiled_seqdep_pitch.data_WithOutlier.(syl_fields{i}){ii,4}; % extract individual specs
    %     end
    %     compiled_seqdep_pitch.summary_stats_WithOutlier.(syl_fields{i}).spec_mean=mean(X,3);
end


%% SAVE

% sequence specific
compiled_seqdep_pitch.PARAMETERS.SEQFILTER.SeqPreList=SeqPreList;
compiled_seqdep_pitch.PARAMETERS.SEQFILTER.SylTargList=SylTargList;
compiled_seqdep_pitch.PARAMETERS.SEQFILTER.SeqPostList=SeqPostList;


% SAVE
timestamp=compiled_seqdep_pitch.PARAMETERS.timestamp; % same as that already used
date=compiled_seqdep_pitch.PARAMETERS.date;
bluejaynum=compiled_seqdep_pitch.PARAMETERS.bluejaynum;
birdname=compiled_seqdep_pitch.PARAMETERS.birdname;
phrase=compiled_seqdep_pitch.PARAMETERS.phrase;

savedir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/compile_seq_dep_pitch_data_' phrase '/SEQFILTER'];
savename=['DataStruct' '_' date];


try
    cd(savedir);
catch err % dir does not exist
    mkdir(savedir);
    cd(savedir);
end

save(savename,'compiled_seqdep_pitch','-v7.3');

cd(curr_dir)

disp('DONE!');

