function spct_strct = make_avg_spects(batch_in,syls,proto_sng,p_f_type,save_name)



for i = 1:length(char(syls))
    syl = char(syls(i));
    [spm,mu_spect,var_spect,axis_data,p_spect,p_snd]=avg_spect_proto_SAM_HACK(char(batch_in),char(syl),char(proto_sng),char(p_f_type));
%    [spm,mu_spect,var_spect,axis_data,p_spect,p_snd]=avg_spect_proto(char(batch_in),char(syl),char(proto_sng),char(p_f_type));
    spct_strct{i}.mu_spct = mu_spect;
    spct_strct{i}.var_spct = var_spect;
    spct_strct{i}.axs_data = axis_data;
    spct_strct{i}.lbl = char(syl);
    spct_strct{i}.proto = p_spect;
    spct_strct{i}.p_wvfrm = p_snd;
    spct_mtrx{i}.spm = spm;

%  [mu_spect, std_spect, axisdata]=avg_spect_KBS_HACK(batch_in);
%     spct_strct{i}.mu_spct = mu_spect;
% %     spct_strct{i}.var_spct = var_spect;
%     spct_strct{i}.axs_data = axis_data;
% %     spct_strct{i}.lbl = char(syl);
% %     spct_strct{i}.proto = p_spect;
% %     spct_strct{i}.p_wvfrm = p_snd;
% %     spct_mtrx{i}.spm = spm;
end

save_name2 = [char(save_name) 'FLOOR.spct_mtrx.mat'];

save(char(save_name2), 'spct_mtrx');

save_name3 = [char(save_name) 'FLOOR.spct_strct.mat'];

save(char(save_name3), 'spct_strct');