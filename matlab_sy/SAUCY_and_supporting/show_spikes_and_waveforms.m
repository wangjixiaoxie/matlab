function foo(ch_num)
n_files=4;
t_per_plot=1;
t_start=2;
show_axis=0;

if isunix
    ls *bin > batchfile
else
    eval(sprintf('!dir /B *.neuralnot_CH%s* > batchfile',num2str(ch_num)))
end

fid=fopen('batchfile','r');
file_ct=1;
spike_color_vec='rgbc';
while file_ct<=n_files
    file_name=fgetl(fid);
    cbin_fname=file_name(1:end-18);
    if (~ischar(file_name))
        fclose(fid);
        break
    end
    disp(['Loading ' file_name]);load(file_name)
    disp(['Loading ' cbin_fname]);[dat,fs]=evsoundin('',cbin_fname,['obs' num2str(ch_num)]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j(file_ct)=subplot(n_files*2,8,(file_ct-1)*16+1);hold on;if ~show_axis;axis off;end
    for x=1:Data.n_units+1;
        ellipsoid=plotcov_sam(Data.cov_centers{x},Data.cov_matrices{x},spike_color_vec(x),.0455);axis tight
    end
    
    h=subplot(n_files*2,8,(file_ct-1)*16+9);hold on;%if ~show_axis;axis off;end
    for x=1:Data.n_units+1;
        plot(Data.mean_waveform{x},'color',spike_color_vec(x),'linew',2);
    end
    
    subplot(n_files,8,((file_ct-1)*8+2):file_ct*8);hold on;if ~show_axis;axis off;end
    if file_ct==1;title(['Plotting ' num2str(t_per_plot) ' sec of waveforms and spikes starting at t=' num2str(t_start)],'fontsize',12,'fontweight','bold');end
    t_lims=[t_start t_start+t_per_plot];
    dat_id=ceil([t_lims(1)*fs:(t_start+t_per_plot)*fs]);
    dat_to_plot=dat(dat_id);
    for x=1:Data.n_units+1;
        spiketimes_cluster=Data.spiketimes{x};
        spiketimes_cluster=spiketimes_cluster(find(spiketimes_cluster>t_lims(1) & spiketimes_cluster<t_lims(2)));
        spiketimes_cluster_id=round((spiketimes_cluster-t_start)*fs)+1;
        plot(spiketimes_cluster_id,dat_to_plot(spiketimes_cluster_id),'ko','markerfacecolor',spike_color_vec(x))
    end
    plot(dat_to_plot,'k')
    axis tight
    xl=get(gca,'xlim');yl=get(gca,'ylim');
    text(mean(xl),yl(2),RemoveUnderScore(cbin_fname),'fontsize',9,'fontweight','bold')
    set(h,'ylim',yl)
    file_ct=file_ct+1;
end

