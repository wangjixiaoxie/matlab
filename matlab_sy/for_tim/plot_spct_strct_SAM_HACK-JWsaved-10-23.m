function plot_spct_strct(spct_strct,dovar)
if nargin == 1, dovar = 0; end

th_low=2.75;th_high=4.5;
use_interpolation=0;

%for i = 3
for i = 1:length(spct_strct)


    SPECT_STRUCT_ADJUSTED=spct_strct{i}.mu_spct;
    id_low=find(SPECT_STRUCT_ADJUSTED<th_low);
    SPECT_STRUCT_ADJUSTED(id_low)=th_low;

    id_high=find(SPECT_STRUCT_ADJUSTED>th_high);
    SPECT_STRUCT_ADJUSTED(id_high)=th_high;

    crntsyl = spct_strct{i}.lbl;
    t = spct_strct{i}.axs_data(1:2); f = spct_strct{i}.axs_data(3:4);
    figure;
    subplot(1,2,1); colormap(hot);
    if use_interpolation
        imagesc(t,f,interp2(SPECT_STRUCT_ADJUSTED,2));
    else
        imagesc(t,f,SPECT_STRUCT_ADJUSTED);
    end
    grid on; title(['MU Spect for ' char(crntsyl)]);
    if use_interpolation
        subplot(1,2,2); imagesc(t,f,interp2(spct_strct{i}.proto,2));
    else
        subplot(1,2,2); imagesc(t,f,spct_strct{i}.proto);
    end
    grid on; title(['Proto Spect for ' char(crntsyl)]);
    if dovar
        figure; imagesc(t,f,spct_strct{i}.var_spct./spct_strct{i}.mu_spct); colormap(hot);colorbar;
    end
    SPECT_STRUCT=spct_strct{i}.mu_spct;
    eval(sprintf('save spect_%s_pre t f SPECT_STRUCT',spct_strct{i}.lbl))

end
