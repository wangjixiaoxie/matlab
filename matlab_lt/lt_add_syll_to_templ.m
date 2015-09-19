%% IN PROGRESS
old_template=input('what is the name of the old template (e.g. pu13bk43b.dat)?', s);
old_template_var=input('retype it, but without .dat')

new_template=input('what is the name of the new template (e.g. pu13bk43b.dat)?', s);

load(['../' old_template])



templ_temporary=[old_template_var 