function [StatsStruct,Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_Statistics(StatsStruct,Params,KeepOutliers);
%% LT 4/7/15 - v2, changed saving format. will also put in wiggle analysis and more complete statistical tests.
% KeepOutliers =1; does not throw out. if 0, then does. 

%% PARAMS
NumTimeWinds=length(Params.TimeField);
TimeWindFields=Params.TimeField;
NumFields=length(Params.FieldsToCheck);
c=1;
alpha=0.05;

%% PERMUTATION TEST TO COMPARE MEANS OF 2 CONDITIONS (FOR ALL MEASURED FEATURES)
disp(' ');
disp('Comparing means using permutation test');

% VARIABLES
Nsims=5000; % number of permutations sims
DataTypeList={'Pitch','Ampl','WEntropy'};

% disp(['Number of resamples: ' num2str(Nsims)]);

for kk=1:length(DataTypeList);
    data_type=DataTypeList{kk};
    
    % RUN
    NumFields_trialtypes=length(Params.FieldsToCheck);
    NumTimeFields=length(Params.TimeField);
    disp(' ');
    disp([data_type ' p values']);
    
    % only perform if am comparing two trial types (majority of cases).
    if NumFields_trialtypes==2;
        
        figure; hold on;
        
        % perform separately for all time windows
        for i=1:NumTimeFields;
            timefield=Params.TimeField{i};
            
            % what datapoints are we comparing?
            if strcmp(data_type,'Ampl')==1;
                X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
                Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
            else
                X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals;
                Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals;
            end
            
            % Remove Outliers?
            if KeepOutliers==0;
                outlinds=StatsStruct.(Params.FieldsToCheck{1}).OutlierInds_AllWindsCombn;
                X(outlinds)=[];
                
                outlinds=StatsStruct.(Params.FieldsToCheck{2}).OutlierInds_AllWindsCombn;
                Y(outlinds)=[];
            end
            
            
            
            % Continue
            Nx=length(X); % sample sizes
            Ny=length(Y);
            
            if c==1;
                disp(' ');
                disp('Sample Sizes:');
                disp([Params.FieldsToCheck{1} ': ' num2str(Nx)]);
                disp([Params.FieldsToCheck{2} ': ' num2str(Ny)]);
                c=2;
            end
            
            D=mean(X)-mean(Y); % real difference
            Z=[X Y]; % combined data
            
            % permute Nsims times
            Dperm=[];
            for ii=1:Nsims;
                
                Iperm=randperm(Nx+Ny); % permute indices
                
                Zperm=Z(Iperm);
                
                Xperm=Zperm(1:Nx);
                Yperm=Zperm(Nx+1:end);
                
                Dperm(ii)=mean(Xperm)-mean(Yperm);
            end
            
            p=length(find(abs(Dperm)>=abs(D)))/Nsims;
            disp([timefield ':  ' num2str(p)]);
            
            
            % PLOT
            subplot(2,NumTimeFields,i); hold on;
            title(['timefield' num2str(i)]);

            Xsem=lt_sem(X);
            Ysem=lt_sem(Y);
            errorbar(i-0.1,mean(X),Xsem,'ko','MarkerFaceColor','r','MarkerSize',7);
            errorbar(i+0.1,mean(Y),Ysem,'ko','MarkerFaceColor','g','MarkerSize',7);
            
             % indicate significance
            Ylims=ylim;
            yrange=Ylims(2)-Ylims(1);
            
            if p<alpha;
                plot(i,max([mean(X)+Xsem mean(Y)+Ysem])+0.015*yrange,'r*','MarkerSize',12);
            end
            
            % plot p-value
            text(i, max([mean(X)+Xsem mean(Y)+Ysem])+0.025*yrange,['p=' num2str(p,'%5.3f')],'Color','b');

            
            % Plot in scatter of datapoints
            subplot(2,NumTimeFields,i+NumTimeFields); hold on;
            
            plot(i-0.12+0.04*rand(length(X),1),X,'or')
            plot(i+0.08+0.04*rand(length(Y),1),Y,'og')
            
            set(gca,'XTick',[1 2]);
            set(gca,'XTickLabel',{Params.FieldsToCheck{1},Params.FieldsToCheck{2}});
            
            
            % % plot histogram
            % figure;
            % hist(Dperm);
            % line([D D],ylim)
            
        end
    end
    lt_subtitle([data_type ' - Means +/- SEM (and permutation test)']);
end


%% PERMUTATION TEST TO COMPARE DIFFERENCE BETWEEN CONTOURS OF 2 CONDITIONS
disp(' ');
disp('Comparing contour differences using permutation test');

% VARIABLES
% Nsims=5000; % number of permutations sims
DataTypeList={'Pitch','Ampl','WEntropy'};

% disp(['Number of resamples: ' num2str(Nsims)]);


for kk=1:length(DataTypeList);
    data_type=DataTypeList{kk};
    if strcmp(data_type,'Pitch')==1;
        % RUN
        NumFields_trialtypes=length(Params.FieldsToCheck);
        NumTimeFields=length(Params.TimeField);
        disp(' ');
        disp([data_type ' p values (differences between mean contour shape of two conditions)']);
        
        % only perform if am comparing two trial types (majority of cases).
        if NumFields_trialtypes==2;
            
            figure; hold on;
            title(['Difference (' Params.FieldsToCheck{2} ' - ' Params.FieldsToCheck{1} ') between mean contour shapes, then mean over time']);
            
            % perform separately for all time windows
            for i=1:NumTimeFields;
                timefield=Params.TimeField{i};
                
                % what samples do I want (of PC)
                ind1_ms=Params.TimeWindowList{i}(1); % in ms
                ind2_ms=Params.TimeWindowList{i}(2);
                
                
                [~,ind1_samp]=min(abs(1000*Params.tf_bins.tPC-ind1_ms)); % convert to sample num
                [~,ind2_samp]=min(abs(1000*Params.tf_bins.tPC-ind2_ms));
                
                
                
                
                % what datapoints are we comparing?
                if strcmp(data_type,'Ampl')==1;
%                     X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
%                     Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
                else
                    
                    
                    X=StatsStruct.(Params.FieldsToCheck{1}).PC(ind1_samp:ind2_samp,:); % matrix of PCs within the window
                    Y=StatsStruct.(Params.FieldsToCheck{2}).PC(ind1_samp:ind2_samp,:);
                end
                
            % Remove Outliers?
            if KeepOutliers==0;
                outlinds=StatsStruct.(Params.FieldsToCheck{1}).OutlierInds_AllWindsCombn;
                X(:,outlinds)=[];
                
                outlinds=StatsStruct.(Params.FieldsToCheck{2}).OutlierInds_AllWindsCombn;
                Y(:,outlinds)=[];
            end
            
            % Continue
                Nx=size(X,2); % sample sizes
                Ny=size(Y,2);
                tmp=abs(mean(X,2)-mean(Y,2)); % take mean across trials, subtract contours and take abs
                D=mean(tmp); % mean abs difference over time
                
                Z=[X Y]; % combined data
                
                % permute Nsims times
                Dperm=[];
                for ii=1:Nsims;
                    
                    Iperm=randperm(Nx+Ny); % permute indices
                    
                    Zperm=Z(:,Iperm);
                    
                    Xperm=Zperm(:,1:Nx);
                    Yperm=Zperm(:,Nx+1:end);
                    
                    tmp=abs(mean(Xperm,2)-mean(Yperm,2)); % subtract contours and take abs
                    Dperm(ii)=mean(tmp); % mean abs difference over time
                    
                end
                
                p=length(find(abs(Dperm)>=abs(D)))/Nsims;
                disp([timefield ':  ' num2str(p)]);
                
                % PLOT
                subplot(2,NumTimeFields,i);
                hold on;
                
                % plot miniature version of contours
                Xstd=std(X,0,2);
                Xsem=Xstd./sqrt(Nx-1);
                Ystd=std(Y,0,2);
                Ysem=Ystd./sqrt(Ny-1);
                
                shadedErrorBar(1:size(X,1),mean(X,2),Xsem,{'Color','r'},1);
                shadedErrorBar(1:size(Y,1),mean(Y,2),Ysem,{'Color','g'},1);
                
                % bar plot of d
                subplot(2,1,2); hold on;
                bar(i, D);
                title('Difference between contours, mean across time');
                ylabel('mean (across time window) contour diff');
                set(gca,'XTick',1:NumTimeFields);
                set(gca,'XTickLabel',Params.TimeField);
                
                % indicate significance
                Ylims=ylim;
                yrange=Ylims(2)-Ylims(1);
                
                if p<alpha;
                    plot(i,D+0.04*yrange,'r*','MarkerSize',12);
                end
                
                % plot p-value
                text(i, D+0.06*yrange,['p=' num2str(p,'%5.3f')],'Color','k');
                
                
                % plot histogram
                %             figure;
                %             hist(Dperm);
                %             line([D D],ylim)
                
            end
            
            lt_subtitle('Pitch contours and difference between contours');
        end
    end
end

%% PERMUTATION TEST TO COMPARE DIFFERENCE BETWEEN CONTOURS (but take diff before averaging)
% mean of one condition, each sample of other condition is compared (i.e.
% each gives a diff value.  take mean of those values. (above is different,
% takes mean across trials first, then compares contours)

if strcmp(Params.FieldsToCheck{2},'StimNotCatch')==0;
    disp('issue, need Params.FieldsToCheck{2} to be StimNotCatch for one thing to work');
    disp('Continuing anyway')
end


disp(' ');
disp('Comparing waveriness differences using permutation test - USING MEAN DIFF');


% VARIABLES
Nsims=5000; % number of permutations sims
DataTypeList={'Pitch','Ampl','WEntropy'};

% disp(['Number of resamples: ' num2str(Nsims)]);


for kk=1:length(DataTypeList);
    data_type=DataTypeList{kk};
    if strcmp(data_type,'Pitch')==1;
        % RUN
        NumFields_trialtypes=length(Params.FieldsToCheck);
        NumTimeFields=length(Params.TimeField);
        disp(' ');
        disp([data_type ' p values']);
        
        figure; hold on;
        
        
        % only perform if am comparing two trial types (majority of cases).
        if NumFields_trialtypes==2;
            
            % perform separately for all time windows
            for i=1:NumTimeFields;
                timefield=Params.TimeField{i};
                
                % what samples do I want (of PC)
                ind1_ms=Params.TimeWindowList{i}(1); % in ms
                ind2_ms=Params.TimeWindowList{i}(2);
                
                
                [~,ind1_samp]=min(abs(1000*Params.tf_bins.tPC-ind1_ms)); % convert to sample num
                [~,ind2_samp]=min(abs(1000*Params.tf_bins.tPC-ind2_ms));
                
                
                
                
                % what datapoints are we comparing?
                if strcmp(data_type,'Ampl')==1;
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
                else
                    X=StatsStruct.(Params.FieldsToCheck{1}).PC(ind1_samp:ind2_samp,:); % matrix of PCs within the window
                    Y=StatsStruct.(Params.FieldsToCheck{2}).PC(ind1_samp:ind2_samp,:);
                end
                
            % Remove Outliers?
            if KeepOutliers==0;
                outlinds=StatsStruct.(Params.FieldsToCheck{1}).OutlierInds_AllWindsCombn;
                X(:,outlinds)=[];
                
                outlinds=StatsStruct.(Params.FieldsToCheck{2}).OutlierInds_AllWindsCombn;
                Y(:,outlinds)=[];
            end
            
            % Continue
                Nx=size(X,2); % sample sizes
                Ny=size(Y,2);
                
                % get actual D
                Xmean=mean(X,2); % mean contour of catch
                
                YminusXmean=Y-repmat(Xmean,1,Ny); % for each Y sample, subtraact mean X contour
                YminusXmean=abs(YminusXmean); % get absolute value
                YminusXmean=mean(YminusXmean,1); % take mean across time for each trial
                
                
                % plot distribution
                %             figure; hold on
                %             title('distribution where each val is one stim trial, controus difference relative to mean non-stim contour');
                %             xlabel('difference FF');
                %             hist(YminusXmean);
                
                % finally take mean across trials
                D=mean(YminusXmean);
                
                % SIMULATE
                Z=[X Y]; % combined data
                
                % permute Nsims times
                Dperm=[];
                for ii=1:Nsims;
                    
                    Iperm=randperm(Nx+Ny); % permute indices
                    
                    Zperm=Z(:,Iperm);
                    
                    Xperm=Zperm(:,1:Nx);
                    Yperm=Zperm(:,Nx+1:end);
                    
                    
                    % Get Dperm
                    Xpermmean=mean(Xperm,2); % mean contour of catch
                    
                    YminusXmeanPerm=Yperm-repmat(Xpermmean,1,Ny); % for each Y sample, subtraact mean X contour
                    YminusXmeanPerm=abs(YminusXmeanPerm); % get absolute value
                    YminusXmeanPerm=mean(YminusXmeanPerm,1); % take mean across time for each trial
                    
                    Dperm(ii)=mean(YminusXmeanPerm);
                    
                end
                
                p=length(find(abs(Dperm)>=abs(D)))/Nsims;
                disp([timefield ':  ' num2str(p)]);
                
                
                % PLOT
                subplot(2,NumTimeFields,i);
                hold on;
                title('Stim contours and mean catch contour');
                
                % plot contours of stim vs. stim catch
                
                plot(Y,'--','Color',[0.85 0.85 0.85]);
                Xstd=std(X,0,2);
                Xsem=Xstd./sqrt(Nx-1);
                shadedErrorBar(1:size(X,1),mean(X,2),Xsem,{'Color','r'},1);
                
                subplot(2,1,2);
                hold on;
                bar(i, D);
                title('Difference between stim contours and mean catch contour');
                ylabel('mean (across time window) contour diff');
                set(gca,'XTick',1:NumTimeFields);
                set(gca,'XTickLabel',Params.TimeField);
                
                % indicate significance
                Ylims=ylim;
                yrange=Ylims(2)-Ylims(1);
                
                if p<alpha;
                    plot(i,D+0.04*yrange,'r*','MarkerSize',12);
                end
                
                % plot p-value
                text(i, D+0.06*yrange,['p=' num2str(p,'%5.3f')],'Color','k');
                
                
                
                %             % plot histogram
                %             figure; hold on;
                %             title([timefield ': mean Difference between each stim trial vs. mean nonstim - Dperm distr. and real D'])
                %             xlabel('integrated FF difference');
                %             hist(Dperm);
                %             line([D D],ylim)
                
            end
            lt_subtitle([data_type ' - Subtract each stim contour from mean catch contour, then mean across time: gives one value']);
            
        end
    end
end

%% DO AGAIN, USING MEDIAN
disp(' ');
disp('Comparing waveriness differences using permutation test - USING MEDIAN DIFF');


% VARIABLES
Nsims=5000; % number of permutations sims
DataTypeList={'Pitch','Ampl','WEntropy'};

% disp(['Number of resamples: ' num2str(Nsims)]);


for kk=1:length(DataTypeList);
    data_type=DataTypeList{kk};
    if strcmp(data_type,'Pitch')==1;
        % RUN
        NumFields_trialtypes=length(Params.FieldsToCheck);
        NumTimeFields=length(Params.TimeField);
        disp(' ');
        disp([data_type ' p values']);
        
        figure; hold on;
        
        
        % only perform if am comparing two trial types (majority of cases).
        if NumFields_trialtypes==2;
            
            % perform separately for all time windows
            for i=1:NumTimeFields;
                timefield=Params.TimeField{i};
                
                % what samples do I want (of PC)
                ind1_ms=Params.TimeWindowList{i}(1); % in ms
                ind2_ms=Params.TimeWindowList{i}(2);
                
                
                [~,ind1_samp]=min(abs(1000*Params.tf_bins.tPC-ind1_ms)); % convert to sample num
                [~,ind2_samp]=min(abs(1000*Params.tf_bins.tPC-ind2_ms));
                
                
                
                
                % what datapoints are we comparing?
                if strcmp(data_type,'Ampl')==1;
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
                else
                    
                    X=StatsStruct.(Params.FieldsToCheck{1}).PC(ind1_samp:ind2_samp,:); % matrix of PCs within the window
                    Y=StatsStruct.(Params.FieldsToCheck{2}).PC(ind1_samp:ind2_samp,:);
                end
    
                
                % Remove Outliers?
                if KeepOutliers==0;
                    outlinds=StatsStruct.(Params.FieldsToCheck{1}).OutlierInds_AllWindsCombn;
                    X(:,outlinds)=[];
                    
                    outlinds=StatsStruct.(Params.FieldsToCheck{2}).OutlierInds_AllWindsCombn;
                    Y(:,outlinds)=[];
                end
            
            % Continue
                Nx=size(X,2); % sample sizes
                Ny=size(Y,2);
                
                % get actual D
                Xmean=mean(X,2); % mean contour of catch
                
                YminusXmean=Y-repmat(Xmean,1,Ny); % for each Y sample, subtraact mean X contour
                YminusXmean=abs(YminusXmean); % get absolute value
                YminusXmean=mean(YminusXmean,1); % take mean across time for each trial
                
                
                % plot distribution
%                 figure; hold on
%                 title('distribution where each val is one stim trial, controus difference relative to mean non-stim contour');
%                 xlabel('difference FF');
%                 hist(YminusXmean);
                
                % finally take mean across trials
                D=median(YminusXmean);
                
                % SIMULATE
                Z=[X Y]; % combined data
                
                % permute Nsims times
                Dperm=[];
                for ii=1:Nsims;
                    
                    Iperm=randperm(Nx+Ny); % permute indices
                    
                    Zperm=Z(:,Iperm);
                    
                    Xperm=Zperm(:,1:Nx);
                    Yperm=Zperm(:,Nx+1:end);
                    
                    
                    % Get Dperm
                    Xpermmean=mean(Xperm,2); % mean contour of catch
                    
                    YminusXmeanPerm=Yperm-repmat(Xpermmean,1,Ny); % for each Y sample, subtraact mean X contour
                    YminusXmeanPerm=abs(YminusXmeanPerm); % get absolute value
                    YminusXmeanPerm=mean(YminusXmeanPerm,1); % take mean across time for each trial
                    
                    Dperm(ii)=median(YminusXmeanPerm);
                    
                end
                
                p=length(find(abs(Dperm)>=abs(D)))/Nsims;
                disp([timefield ':  ' num2str(p)]);
                
                % PLOT
                subplot(2,NumTimeFields,i);
                hold on;
                title('Stim contours and mean catch contour');
                
                % plot contours of stim vs. stim catch
                
                plot(Y,'--','Color',[0.85 0.85 0.85]);
                Xstd=std(X,0,2);
                Xsem=Xstd./sqrt(Nx-1);
                shadedErrorBar(1:size(X,1),mean(X,2),Xsem,{'Color','r'},1);
                
                subplot(2,1,2);
                hold on;
                bar(i, D);
                title('Difference between stim contours and mean catch contour');
                ylabel('MEDIAN diff between stim contours and mean catch contour');
                set(gca,'XTick',1:NumTimeFields);
                set(gca,'XTickLabel',Params.TimeField);
                
                % indicate significance
                Ylims=ylim;
                yrange=Ylims(2)-Ylims(1);
                
                if p<alpha;
                    plot(i,D+0.04*yrange,'r*','MarkerSize',12);
                end
                
                % plot p-value
                text(i, D+0.06*yrange,['p=' num2str(p,'%5.3f')],'Color','k');
                
                % plot histogram
%                 figure; hold on;
%                 title([timefield ': mean Difference between each stim trial vs. mean nonstim - Dperm distr. and real D'])
%                 xlabel('integrated FF difference');
%                 hist(Dperm);
%                 line([D D],ylim)
                
            end
        end
    end
end
%% COMPARE WAVERINESS - IN PROGRESS
% each condition subtract its own mean.
% IN PROGRESS - should subtract each trial's mean (not global mean).

if (0)
disp(' ');
disp('Comparing waveriness differences (from itself) using permutation test - USING MEAN DIFF');


% VARIABLES
Nsims=5000; % number of permutations sims
DataTypeList={'Pitch','Ampl','WEntropy'};

% disp(['Number of resamples: ' num2str(Nsims)]);


for kk=1:length(DataTypeList);
    data_type=DataTypeList{kk};
    if strcmp(data_type,'Pitch')==1;
        % RUN
        NumFields_trialtypes=length(Params.FieldsToCheck);
        NumTimeFields=length(Params.TimeField);
        disp(' ');
        disp([data_type ' p values (waveriness vs. own mean)']);
        
        figure; hold on;
        title([data_type ' - waveriness vs. own mean']);
        
        
        % only perform if am comparing two trial types (majority of cases).
        if NumFields_trialtypes==2;
            
            % perform separately for all time windows
            for i=1:NumTimeFields;
                timefield=Params.TimeField{i};
                
                % what samples do I want (of PC)
                ind1_ms=Params.TimeWindowList{i}(1); % in ms
                ind2_ms=Params.TimeWindowList{i}(2);
                
                
                [~,ind1_samp]=min(abs(1000*Params.tf_bins.tPC-ind1_ms)); % convert to sample num
                [~,ind2_samp]=min(abs(1000*Params.tf_bins.tPC-ind2_ms));
                
                
                
                
                % what datapoints are we comparing?
                if strcmp(data_type,'Ampl')==1;
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
                else
                    X=StatsStruct.(Params.FieldsToCheck{1}).PC(ind1_samp:ind2_samp,:); % matrix of PCs within the window
                    Y=StatsStruct.(Params.FieldsToCheck{2}).PC(ind1_samp:ind2_samp,:);
                end
                
                Nx=size(X,2); % sample sizes
                Ny=size(Y,2);
                
                % get actual D
                
                Xmean=mean(X,2); % mean contour of catch
                Ymean=mean(Y,2);
                
                
                YMinusMean=Y-repmat(Ymean,1,Ny); % for each Y sample, subtraact mean Y contour
                YMinusMeanAbs=abs(YMinusMean); % get absolute value
                YMinusMeanAbs_mean=mean(YMinusMeanAbs,1); % take mean across time for each trial
                
                XMinusMean=X-repmat(Xmean,1,Nx); % for each Y sample, subtraact mean Y contour
                XMinusMeanAbs=abs(XMinusMean); % get absolute value
                XMinusMeanAbs_mean=mean(XMinusMeanAbs,1); % take mean across time for each trial
                
                
                % finally take mean across trials
                D=mean(YMinusMeanAbs_mean)-mean(XMinusMeanAbs_mean);
                
                % SIMULATE
                Z=[X Y]; % combined data
                
                % permute Nsims times
                Dperm=[];
                for ii=1:Nsims;
                    
                    Iperm=randperm(Nx+Ny); % permute indices
                    
                    Zperm=Z(:,Iperm);
                    
                    Xperm=Zperm(:,1:Nx);
                    Yperm=Zperm(:,Nx+1:end);
                    
                    
                    % Get Dperm
                    Xmeanperm=mean(Xperm,2); % mean contour of catch
                    Ymeanperm=mean(Yperm,2);
                    
                    
                    YMinusMeanperm=Yperm-repmat(Ymeanperm,1,Ny); % for each Y sample, subtraact mean Y contour
                    YMinusMeanAbsperm=abs(YMinusMeanperm); % get absolute value
                    YMinusMeanAbs_meanperm=mean(YMinusMeanAbsperm,1); % take mean across time for each trial
                    
                    XMinusMeanperm=Xperm-repmat(Xmeanperm,1,Nx); % for each Y sample, subtraact mean Y contour
                    XMinusMeanAbsperm=abs(XMinusMeanperm); % get absolute value
                    XMinusMeanAbs_meanperm=mean(XMinusMeanAbsperm,1); % take mean across time for each trial
                    
                    
                    % finally take mean across trials
                    Dperm(ii)=mean(YMinusMeanAbs_meanperm)-mean(XMinusMeanAbs_meanperm);
                    
                end
                
                p=length(find(abs(Dperm)>=abs(D)))/Nsims;
                disp([timefield ':  ' num2str(p)]);
                
                
                % PLOT
                
                % plot contours of stim vs. stim catch
                subplot(3,NumTimeFields,i);
                hold on; title('Stim');
                
                
                plot(YMinusMean+mean(Ymean),'-','Color',[0.85 0.85 0.85]);
                plot(Ymean,'-','Color','g');
                
                % plot catch
                subplot(3,NumTimeFields,i+NumTimeFields);
                hold on; title('Catch');
                
                plot(XMinusMean+mean(Xmean),'-','Color',[0.85 0.85 0.85]);
                plot(Xmean,'-','Color','r');
                
                
                
                
                
                
                %             % plot histogram
                %             figure; hold on;
                %             title([timefield ': mean Difference between each stim trial vs. mean nonstim - Dperm distr. and real D'])
                %             xlabel('integrated FF difference');
                %             hist(Dperm);
                %             line([D D],ylim)
                
            end
        end
    end
end

end


%% PERMUTATION TEST TO COMPARE MEDIANS OF 2 CONDITIONS (FOR ALL MEASURED FEATURES)
% outliers not removed
if (0)
    disp(' ');
    disp('Comparing medians using permutation test');
    
    % VARIABLES
    Nsims=5000; % number of permutations sims
    DataTypeList={'Pitch','Ampl','WEntropy'};
    
    disp(['Number of resamples: ' num2str(Nsims)]);
    
    for kk=1:length(DataTypeList);
        data_type=DataTypeList{kk};
        
        % RUN
        NumFields_trialtypes=length(Params.FieldsToCheck);
        NumTimeFields=length(Params.TimeField);
        disp(' ');
        disp([data_type ' p values']);
        
        % only perform if am comparing two trial types (majority of cases).
        if NumFields_trialtypes==2;
            
            % perform separately for all time windows
            for i=1:NumTimeFields;
                timefield=Params.TimeField{i};
                
                % what datapoints are we comparing?
                if strcmp(data_type,'Ampl')==1;
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
                else
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals;
                end
                
                Nx=length(X); % sample sizes
                Ny=length(Y);
                D=median(X)-median(Y); % real difference
                Z=[X Y]; % combined data
                
                % permute Nsims times
                Dperm=[];
                for ii=1:Nsims;
                    
                    Iperm=randperm(Nx+Ny); % permute indices
                    
                    Zperm=Z(Iperm);
                    
                    Xperm=Zperm(1:Nx);
                    Yperm=Zperm(Nx+1:end);
                    
                    Dperm(ii)=median(Xperm)-median(Yperm);
                end
                
                p=length(find(abs(Dperm)>=abs(D)))/Nsims;
                disp([timefield ':  ' num2str(p)]);
                
                % % plot histogram
                % figure;
                % hist(Dperm);
                % line([D D],ylim)
                
            end
        end
    end
    
end
%% PERMUTATION TEST TO COMPARE MEANS OF 2 CONDITIONS (FOR ALL MEASURED FEATURES)
% outliers not yet removed

if (0)
    disp(' ');
    disp('Comparing means using permutation test');
    
    % VARIABLES
    Nsims=5000; % number of permutations sims
    DataTypeList={'Pitch','Ampl','WEntropy'};
    
    disp(['Number of resamples: ' num2str(Nsims)]);
    
    for kk=1:length(DataTypeList);
        data_type=DataTypeList{kk};
        
        % RUN
        NumFields_trialtypes=length(Params.FieldsToCheck);
        NumTimeFields=length(Params.TimeField);
        disp(' ');
        disp([data_type ' p values']);
        
        % only perform if am comparing two trial types (majority of cases).
        if NumFields_trialtypes==2;
            
            % perform separately for all time windows
            for i=1:NumTimeFields;
                timefield=Params.TimeField{i};
                
                % what datapoints are we comparing?
                if strcmp(data_type,'Ampl')==1;
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals_log;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals_log;
                else
                    X=StatsStruct.(Params.FieldsToCheck{1}).WINDOWED.(timefield).(data_type).vals;
                    Y=StatsStruct.(Params.FieldsToCheck{2}).WINDOWED.(timefield).(data_type).vals;
                end
                
                Nx=length(X); % sample sizes
                Ny=length(Y);
                D=mean(X)-mean(Y); % real difference
                Z=[X Y]; % combined data
                
                % permute Nsims times
                Dperm=[];
                for ii=1:Nsims;
                    
                    Iperm=randperm(Nx+Ny); % permute indices
                    
                    Zperm=Z(Iperm);
                    
                    Xperm=Zperm(1:Nx);
                    Yperm=Zperm(Nx+1:end);
                    
                    Dperm(ii)=mean(Xperm)-mean(Yperm);
                end
                
                p=length(find(abs(Dperm)>=abs(D)))/Nsims;
                disp([timefield ':  ' num2str(p)]);
                
                
                % % plot histogram
                % figure;
                % hist(Dperm);
                % line([D D],ylim)
                
            end
        end
    end
    
end

%% BOOTSTRAP TO GET ERROR BARS FOR STD + PLOT
disp('Plotting STDs');
disp('doing bootstrap to get confidence intervals for std');

Nsamps=5000; % number of times to perform sampling
Qfunc='std';

figure; hold on;
PlotColors=lt_make_plot_colors(NumFields,0);

for i=1:NumTimeWinds;
    timewind=TimeWindFields{i};
    subplot(1,NumTimeWinds,i); hold on;
    title(timewind); ylabel('Hz');
    
    for ii=1:NumFields; % for each trail type
        fieldname=Params.FieldsToCheck{ii};

        X=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.vals; % data points
        
        % Remove Outliers?
        if KeepOutliers==0;
            outlinds=StatsStruct.(fieldname).OutlierInds_AllWindsCombn;
            X(outlinds)=[];
        end
        
        Xstd=std(X);
        
        % Bootstrap
        BootStats=lt_bootstrap(X,Qfunc,Nsamps);
        StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.BootStats.STD=BootStats;
        
        % PLOT
        hfig5(ii)=plot(0.5*ii-0.25,Xstd,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
        
        % plot 95% ci
        y1=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.BootStats.STD.CI(1);
        y2=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.BootStats.STD.CI(2);
        
        line([0.5*ii-0.25 0.5*ii-0.25],[y1 y2]);
    end
end

legend(hfig5,Params.FieldsToCheck);
lt_subtitle('Pitch STD (95% CI), sorted by time window and trial');


%% PLOT STD - I put into above
% PITCH - median
% figure; hold on;
% PlotColors=lt_make_plot_colors(NumFields,0);
% for i=1:NumTimeWinds;
%     timewind=TimeWindFields{i};
%     subplot(1,NumTimeWinds,i); hold on;
%     title(timewind); ylabel('Hz');
%     
%     for ii=1:NumFields; % for each trail type
%         fieldname=Params.FieldsToCheck{ii};
%         
%         % plot medain vals
%         Y=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SD;
%         
%         % Remove Outliers?
%         if KeepOutliers==0;
%             outlinds=StatsStruct.(fieldname).OutlierInds_AllWindsCombn;
%             X(outlinds)=[];
%         end
%         
%         %         hfig5(ii)=plot(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.BootStats.STD.MEDIAN,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
%         hfig5(ii)=plot(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SD,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
%         
%         % plot 95% conf intervals
%         y1=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.BootStats.STD.CI(1);
%         y2=StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.BootStats.STD.CI(2);
%         
%         line([0.5*ii-0.25 0.5*ii-0.25],[y1 y2]);
%     end
% end
% 
% legend(hfig5,Params.FieldsToCheck);
% lt_subtitle('Pitch STD (95% CI), sorted by time window and trial');
% 

%% COMPARE VARIANCE



% Pitch wiggle



% variability across trials
% PITCH - median
% figure; hold on;
% for i=1:NumTimeWinds;
%     timewind=TimeWindFields{i};
%     subplot(1,NumTimeWinds,i); hold on;
%     title(timewind); ylabel('Hz');
%
%     for ii=1:NumFields; % for each trail type
%         fieldname=Params.FieldsToCheck{ii};
%
%         % plot mean vals
%         hfig5(ii)=plot(0.5*ii-0.25,StatsStruct.(fieldname).WINDOWED.(timewind).Pitch.SD,'s','Color','k','MarkerFaceColor',PlotColors{ii}, 'MarkerSize',8);
%     end
% end
%
% legend(hfig5,Params.FieldsToCheck);
% lt_subtitle('Pitch STD (MEDIAN/SEM), sorted by time window and trial');
%



% Broken correlations acro



%% SAVE

disp('Saving...')
cd(Params.savefolder);
tstamp=lt_get_timestamp(0);

% done note
DoneNote=['DONE_Statistics_' tstamp '.txt'];
fid1=fopen(DoneNote,'w');
fclose(fid1);

% save figures
try cd('FIGURES/Statistics');
catch err
mkdir('FIGURES/Statistics');
cd('FIGURES/Statistics');
end

lt_save_all_figs
cd('../../');

disp('Done!');





