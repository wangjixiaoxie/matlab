% %getstyle2.m  9.20.08
% %used to set style for inactivplotcompar

%possible style inputs are 'co' 'do' 'up'


% % 1. ps.color
% % 2. ps.width
% % 3. ps.marksize
% % 4. ps.markcolor
% % 5. ps.markertype
% 
% color, width of line,
% %marker size, color of line
% 

function [ps,ms] = getstyle2(style)



if (style=='do') 
    ps(1).col='r'
    ps(1).wid=2;
    ps(1).fillcol{1}=[0.82 0.82 0.82]
    ps(1).fillcol{2}=[0.92 0.96 0.98]
    ms=[];;
   
    
elseif(style=='up')
    ps(1).col='k'
    ps(1).fillcol{1}=[0.82 0.82 0.82]
    ps(1).fillcol{2}=[0.92 0.96 0.98]
    ps(1).wid=2;
    ms=[];
elseif(style=='co')
    ps(1).col='k'
    ps(1).wid=2;
    ps(1).fillcol{1}=[0.82 0.82 0.82]
    ps(1).fillcol{2}=[0.92 0.96 0.98]
    ms=[];
elseif(style=='ex')
    ps(1).col=[0 0 0]
    ps(1).wid=2;
    ps(2).col=[0.4 0.4 1]
    ps(2).wid=2;
    ps(1).fillcol=[0.82 0.82 0.82]
    ps(2).fillcol=[0.92 0.96 0.98]
end
    
%marker1 is black
%marker2 is red.

%for control plots
%line is green
%markersize 1 is black
%markersize 2 is red