function h = imagesc1(x,y,a,clim)
%IMAGESC1 Scale data and display as image.
%       IMAGESC(...) is the same as IMAGE(...) except the data is scaled
%       to use the full colormap.
%       An optional final argument CLIMS = [CLOW CHIGH] can specify the
%       scaling.
%
%       Note: this function ins the same as matlabs imagesc except that it has been modified to work
%             correctly in the case where 'a' has negative entries and CLIMS is arbitrary 
%             (formerly screwed this up, producing 0 or negative indices into colormap)
%              
%             Also have disabled the use of userdata property so that COLORBAR can no longer annotate 
%             meaningful labels.
%
%

%       to re-enable this feature, restore commented out line 
%       NOTE:  IMAGESC places scaling information in the displayed image's
%              UserData property so the COLORBAR command can annotate
%              meaningful labels .
%
%       See also IMAGE, COLORBAR.

%       Copyright (c) 1984-94 by The MathWorks, Inc.
 


cm = colormap;
m = size(cm,1);
if nargin <= 2
        a = x;
end
if rem(nargin,2)  %case where autoscale to data min and max
   amin = min(min(a));
   amax = max(max(a));   
elseif nargin == 4  %case where all arguments supplied
   amin = clim(1);
   amax = clim(2);   
elseif nargin == 2  %case where only data and clim supplied 
   amin = y(1);
   amax = y(2);   
end

idx=max(1,ceil(m*(max(amin,min(amax,a))-amin)/(amax-amin)));

if nargin <= 2
        hh = image(idx);
else
        hh = image(x,y,idx);
end
colormap(cm);
%set(hh,'UserData',[amin amax]);
if nargout > 0
        h = hh;
end
