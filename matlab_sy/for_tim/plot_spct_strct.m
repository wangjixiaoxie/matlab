function plot_spct_strct(spct_strct,dovar)
if nargin == 1, dovar = 0; end

for i = 1:length(spct_strct)
    crntsyl = spct_strct{i}.lbl; 
    t = spct_strct{i}.axs_data(1:2); f = spct_strct{i}.axs_data(3:4); 
    figure; subplot(1,2,1); colormap(hot); 
    imagesc(t,f,interp2(spct_strct{i}.mu_spct,2));
    grid on; title(['MU Spect for ' char(crntsyl)]);
    subplot(1,2,2); imagesc(t,f,interp2(spct_strct{i}.proto,2));
    grid on; title(['Proto Spect for ' char(crntsyl)]);
    if dovar
        figure; imagesc(t,f,spct_strct{i}.var_spct./spct_strct{i}.mu_spct); colormap(hot);colorbar;
    end
end

