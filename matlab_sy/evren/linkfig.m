function linkfig(figid,linkstr);
%linkfig(figid,linkstr);
%link up all the subplots on a given figure

if (exist('linkstr'))
	if (length(linkstr)>0)
		linkaxes(get(figid,'Children'),linkstr);
	else
		linkaxes(get(figid,'Children'));
	end
else
	linkaxes(get(figid,'Children'));
end
return;
