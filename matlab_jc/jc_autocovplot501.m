function jc_autocovplot501(arrayD,arrayUD,arrayL,first_note,last_note,start_bin,stop_bin)
figure; hold on
subplot(131); hold on
jc_autocov501(arrayD,first_note,last_note,start_bin,stop_bin)
subplot(132); hold on
jc_autocov501(arrayUD,first_note,last_note,start_bin,stop_bin)
subplot(133); hold on
jc_autocov501(arrayL,first_note,last_note,start_bin,stop_bin)