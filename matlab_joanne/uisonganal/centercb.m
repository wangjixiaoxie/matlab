function centercb(action,varargin)

% Need button up and down functions here
% disp(action);

switch action
  case 'center_pointer'
    % disp('Got to center_pointer');
    % oldptr = get(gcbf,'Pointer');
    POINTER_IN_AXES = 0;
    fig = get(0,'PointerWindow'); 
    % Look for quick exit
    if fig==0,
      return
    end    
    
    saveunits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    figpos = get(gcbf,'Position');
    ha = findobj(gcbf,'Type','axes','Tag','SAAxis');
    for i=1:length(ha)
      ptrloc= get(0,'PointerLocation');
      ptrloc = (ptrloc - figpos(1,1:2))./figpos(1,3:4);
      % If pointer is in figure then check if it is in an axis
      if (ptrloc(1) <= 1 & ptrloc(1) >= 0 & ptrloc(2) <= 1 & ptrloc(2) >= 0)
	ax = ha(i);
	saveaxesunits = get(ax,'Units');
	set(ax,'Units','normalized');
	axispos = get(ax,'Position');
	if (axispos(1) <= ptrloc(1) & ptrloc(1) <= axispos(1)+axispos(3) & ...
	      axispos(2) <= ptrloc(2) & ptrloc(2) <= axispos(2)+axispos(4))	  
	  % disp('Got to pointer change');
	  POINTER_IN_AXES = 1;
	  set(ax,'Units',saveaxesunits);
	  break;
	else
	  set(ax,'Units',saveaxesunits);
	end % if	
      end
    end

    if POINTER_IN_AXES
      set(gcbf,'Pointer','cross');
%      setptr(gcbf,'glass');
    else
      set(gcbf,'Pointer','arrow');
    end
    set(gcbf,'Units',saveunits);
    
  case 'center'
    POINTER_IN_AXES = 0;
    fig = get(0,'PointerWindow'); 
    % Look for quick exit
    if fig==0,
      return
    end    
    
    saveunits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    figpos = get(gcbf,'Position');
    ha = findobj(gcbf,'Type','axes','Tag','SAAxis');
    for i=1:length(ha)
      ptrloc= get(0,'PointerLocation');
      ptrloc = (ptrloc - figpos(1,1:2))./figpos(1,3:4);
      % If pointer is in figure then check if it is in an axis
      if (ptrloc(1) <= 1 & ptrloc(1) >= 0 & ptrloc(2) <= 1 & ptrloc(2) >= 0)
	ax = ha(i);
	saveaxesunits = get(ax,'Units');
	set(ax,'Units','normalized');
	axispos = get(ax,'Position');
	if (axispos(1) <= ptrloc(1) & ptrloc(1) <= axispos(1)+axispos(3) & ...
	      axispos(2) <= ptrloc(2) & ptrloc(2) <= axispos(2)+axispos(4))	  
	  % disp('Got to pointer change');
	  POINTER_IN_AXES = 1;
	  set(ax,'Units',saveaxesunits);
	  break;
	else
	  set(ax,'Units',saveaxesunits);
	end % if	
      end
    end

    if POINTER_IN_AXES
      pt = get_currentpoint(ax);
      xlim = get_xlim(ax);    
      % Center axis at pt if possible
      xrange = xlim(2) - xlim(1);
      xshift = pt(1) - (xlim(1) + xrange/2);
      xlim_new = xlim + xshift;
      % Fix this up for log scale
      for i=1:length(ha)
	set(ha(i),'XLim',xlim_new);
      end      
    end

    set(gcbf,'Units',saveunits);

    
    
    
    
    
%    % Maybe have some settings to prevent flashing 
%    % when zoom hasn't changed. Must keep prev state.
%    fac = varargin{1};
%    val = get(gcbo,'Value');
%    % Unset all other radiobuttons in group
%    hrb = findobj(gcbf,'Tag','ZoomRb');
%    zrbinfo = str2num(char(get(hrb,'String')));
%    set(hrb(zrbinfo~=fac),'Value',0);
%    ha = findobj(gcbf,'Type','axes');
%    for i=1:length(ha)
%      typ = get(get(ha(i),'Children'),'Type');
%      if strcmp(typ,'image') 
%	axes(ha(i));
%	myzoom out;
%	if val
%	  myzoom(fac);
%	end
%      end
%    end
    
%  case 'zoom_mouse'
%    % disp('Got to zoom_mouse');
%    val = get(gcbo,'Value');
%    if val == 1
%      set(gcbf,'WindowButtonMotionFcn','zoomcb(''zoom_pointer'')');
%      myzoom(gcbf,'on');
%      % Set pointer to mag in all axes
%      % Hold zoom button down, disable all others
%      PushUpOtherTB(gcbo);
%    else
%      zoomcb('zoom_off');
%    end
%  case 'zoom_off'
%    % disp('Got here');
%    myzoom(gcbf,'off');
%    % Reset pointer
%    set(gcbf,'Pointer','arrow');
%    set(gcbf,'WindowButtonMotionFcn','');
%    % Push zoom button up if pressed    
%    htb = findobj(gcbf,'Tag','ZoomTb');
%    set(htb,'Value',0);
%    
%  case 'zoom_out'
%    ha = findobj(gcbf,'Type','axes');
%    for i=1:length(ha)
%      typ = get(get(ha(i),'Children'),'Type');
%      if strcmp(typ,'image') == 1 
%	axes(ha(i));
%	myzoom out;
%      end
%    end
end    
    

function p = get_currentpoint(ax)
%GET_CURRENTPOINT Return equivalent linear scale current point
p = get(ax,'currentpoint'); p = p(1,1:2);
if strcmp(get(ax,'XScale'),'log'),
   p(1) = log10(p(1));
end
if strcmp(get(ax,'YScale'),'log'),
   p(2) = log10(p(2));
end

function xlim = get_xlim(ax)
%GET_XLIM Return equivalent linear scale xlim
xlim = get(ax,'xlim');
if strcmp(get(ax,'XScale'),'log'),
   xlim = log10(xlim);
end

function ylim = get_ylim(ax)
%GET_YLIM Return equivalent linear scale ylim
ylim = get(ax,'ylim');
if strcmp(get(ax,'YScale'),'log'),
   ylim = log10(ylim);
end
