function b = isint(x)
%ISINT    True for integers.
%   ISINT(X) returns 1's where the elements of X are integers and 0's
%   where they are not. For example, ISINT([ 2e9 pi 3+5i NaN ])
%   is [ 1 0 1 0 ].

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-10-17 14:01:15
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(1, 1, nargin));
   if ~isnumeric(x)
      error('Argument must be a numeric array.');
   end

   cls = class(x);                      % class of input argument
   if isempty(x)
      b = feval(cls, x);                % return empty array of same class
   else
      switch cls
         case 'double'
            b = x == round(x);
         case 'single'
            % convert input to double, compare, and convert back
            d = double(x);
            b = single(d == round(d));
         case {'int8', 'int16', 'int32', 'uint8', 'uint16', 'uint32'}
            % return an array of ones of the same class as input argument
            b = repmat(feval(cls, 1), size(b));
         otherwise
            error('Argument is of unrecognized class.');
      end
   end
