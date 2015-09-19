function [elts,nelt]=elements(x)
%ELEMENTS Elements of a matrix in set theoretic sense.
% ELEMENTS(X) is a row vector containing the
% distinct elements of X.
% [ELTS,NELT] = ELEMENTS(X) produces
% the elements of X in ELTS and
% the corresponding count of the elements in NELT.
% X is treated as 1 set and may contain NaN's and Inf's
% (which are counted also). Complex arrays as well as
% sparse matrices and text strings are handled properly.
% ELTS is sorted.
% Enter 'elements' for a demo.

% Author:  J. Rodney Jee, rodjee@delphi.com,  28-JAN-95

if isstr(x)
   str_flag=1;
else str_flag=0;
end 
    
    
if ( nargin == 0 )               % DEMO this function when no input.
   disp('+++++++++++++++ DEMO of ELEMENTS +++++++++++++++')
   disp('GIVEN a set, say')
   x = round( rand(4,6)*4 )
   disp('ELEMENTS returns its')
   [members,counts]=elements(x);
   members
   disp('and their respective')
   counts
   disp('+++++++++++++++++ END of DEMO +++++++++++++++++')
   return
end


% The key ideas of this method are to (1)sort the data, (2)take the
% differences of the sorted data and look for nonzeros in the 
% differences which mark the ends of strings of the same values, (3)
% collect the values of step 2, and (4)use the indices of the
% jumps to tally the members.

if (issparse(x))                         % Check for sparse matrix.
   nzeros=prod(size(x))-nnz(x);
   x = nonzeros(x);                      % Required to be a column matrix.
else
   if ( isstr(x) )                       % Convert text strings to integer.
      xstring=1;
      x = abs(x);
	  x = x(:);
      nzeros=1;                 %added without deep understanding 4/3/03 to fix matlab6.5 failure
   else
      xstring=0;
      nzeros=0;
      x = x(:);
   end
end

indexf = finite(x);
xout   = x( ~indexf );                   % Set aside NaNs and Infs.
x      = sort( x(indexf) );                             % Step (1).
if ( isempty(x) )
   elts = [];
   nelt = [];
elseif (length(x) == 1)
   elts = x;
   nelt = 1;
else
  indjump = find(diff(x) ~= 0);                         % Step (2).
  if ( length(indjump) == 0 )
     elts = x(1);
     nelt = length(x);
  else
     elts = x(indjump);                                 % Step (3).
     elts = [elts.',x(length(x))];
     nelt = diff( [0, indjump.', length(x)] );          % Step (4).
  end
end

if (isempty(xout) & (nzeros==0))        
   if str_flag==1
       elts=setstr(elts);
   end        
    return
else                                  % Append NaN's,Inf's, and 0's.
   nnan = sum(isnan(xout));
   ninf = length(xout) - nnan;
   if ( nnan > 0)
      elts = [elts, NaN];
      nelt = [nelt, nnan];
   end
   if ( ninf > 0)
      elts = [elts, Inf];
      nelt = [nelt, ninf];
   end
   if ( nzeros > 0)
      elts = [elts, 0];
	  nelt = [nelt, nzeros];
   end
end

if (xstring)                        
   elts = setstr(elts);
end

if str_flag==1
   elts=setstr(elts);
end   
