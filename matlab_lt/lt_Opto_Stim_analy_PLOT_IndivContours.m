%% LT 2/21/15 - plots individual renditions, so you can look for outliers, etc.
% TO DO: 
% 1) make title tell me what song file it is
% 2) and/or plot spectrogram 
% 

function lt_Opto_Stim_analy_PLOT_IndivContours(Params, StatsStruct,fieldname,bins);



fieldname='COMBO_StimNotCatch_TargAllCatch';
bins=[220 300];

%% PARAMS
NumFields=length(Params.FieldsToCheck);
subs_per_fig=16;


%% PLOT PITCH CONTOURS IN ORDER

tPC=Params.tf_bins.tPC;
tSP=Params.tf_bins.tSP;
fSP=Params.tf_bins.fSP;


t1=find(tPC-bins(1)/1000>0,1,'first'); % start index
t2=find(tPC-bins(2)/1000<0,1,'last'); % last index
    
%     fieldname=Params.FieldsToCheck{iii};
    NumSyls=size(StatsStruct.(fieldname).PC,2);
    
    [num_figures, row_plots, col_plots]=lt_get_subplot_size(NumSyls,subs_per_fig); % how many plots?
    hold on;
    PC=StatsStruct.(fieldname).PC;
    c=1; % initiate counter
    
    for j=1:NumSyls;
        
    fignum=ceil((c-0.5)/subs_per_fig);
    figure(fignum); hold on;

    % 1) PC
    
    if mod(c,subs_per_fig)==0;
        splotNum=subs_per_fig;
    else
        splotNum=mod(c,subs_per_fig);
    end
    
    
    subplot(4,4,splotNum); hold on;
    plot(tPC(t1:t2)*1000,PC(t1:t2,j),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
    plot(tPC(t1:t2)*1000,mean(PC((t1:t2),:)'),'k')
    ylabel('Frequency (Hz)')
%     title([  );
    xlim([tPC(t1)*1000 tPC(t2)*1000])
    
    c=c+1;
end
%     % 2) SPEC
%     sp_mean=StatsStruct.(fieldname).sp_mean;
%     
%     
%     % first, convert any sp values of 0 to non-zero(to the lowest value present);
%     % solves problem of taking log of 0
%     pp=find(sp_mean>0);
%     mntmp = min(min(sp_mean(pp)));
%     pp=find(sp_mean==0);
%     sp_mean(pp) = mntmp;
%     
%     % second, take log
%     sptemp=log(sp_mean);
%     sptemp = sptemp - min(min(sptemp));
%     sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
%     
%     h1(2)=subplot(4,1,2); hold on;
%     imagesc(tSP*1000, fSP, sptemp);
%     ylabel('Frequency (hz)');
%     axis([tSP(1) tSP(end) fSP(1) fSP(end)]);
%     
%     
end