%plterrbar.m
%This program plots bars with error bar...in stimf format
% filename_arr{1}='2672';
% filename_arr{2}='2672m2';
% filename_arr{3}='2672rev';
% filename_arr{4}='2672m12';
% filename_arr{5}='2672p12';
% filename_arr{6}='2672sn';
clusts_to_plot=[1 ]
stim_to_plt=[1:6]

for ii=1:length(clusts_to_plot)
    subplot(length(clusts_to_plot),1,ii);
    meanspikerat=[];
    errspikerat=[];
    for jj=1:length(stim_to_plt)
        meanspikerat=[meanspikerat stimf(ii,jj).meansongspks];
        errspikerat=[errspikerat stimf(ii,jj).stderrsongspks];
        
        errorbar(meanspikerat,errspikerat,'+','Linewidth',3)
        set(gca,'xtick',[1:length(stim_to_plt)]);
        box off;
        if(ii~=length(clusts_to_plot))
            set(gca,'xcolor','w')
        end
        %set(gca,'xticklabel',['bos';'rev';'l  ';'m10';'m10';'m5 ';'p5 ';'m1t';'p1t';'qui']);
        if(ii==length(clusts_to_plot))
            set(gca,'xticklabel',['bos';'m2 ';'rev';'m12';'p12'; 'sn '],'Fontsize',16);
        end
        YLabel('avg song spike rate (Hz)','Fontsize',16) ;
    end
end
    