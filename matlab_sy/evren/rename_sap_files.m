function rename_sap_files(ff,bird_name,year);
%rename_sap_files(ff,bird_name,year);
%

for ii=1:length(ff)
	fn=ff(ii).name;

	if (~exist(fn,'file'))
		continue;
	end

        tmp=findstr(fn,'_');
        mon=fn(tmp(2)+1:tmp(3)-1);
        if length(mon)<2; mon=['0',mon];end

        dy=fn(tmp(3)+1:tmp(4)-1);
        if length(dy)<2; dy=['0',dy];end

        hr=fn(tmp(4)+1:tmp(5)-1);
        if length(hr)<2; hr=['0',hr];end

        mn=fn(tmp(5)+1:tmp(6)-1);
        if length(mn)<2; mn=['0',mn];end

        tmp2=findstr(fn,'.wav');
        sc=fn(tmp(6)+1:tmp2(1)-1);
        if length(sc)<2; sc=['0',sc];end

        %get an evtaf style name
        fn_new=[bird_name,'_',dy,mon,year,'_',hr,mn,sc,'.wav'];

	eval(['!mv ',fn,' ',fn_new]);
	%disp(['!mv ',fn,' ',fn_new]);
	
end
return
