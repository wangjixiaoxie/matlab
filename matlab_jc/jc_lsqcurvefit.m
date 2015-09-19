 function [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqcurvefit(FUN,x,XDATA,YDATA,LB,UB,options,varargin)
%LSQCURVEFIT solves non-linear least squares problems.
%   LSQCURVEFIT attempts to solve problems of the form:
%   min  sum {(FUN(X,XDATA)-YDATA).^2}  where X, XDATA, YDATA and the values 
%    X                                  returned by FUN can be vectors or 
%                                       matrices.
%
%   X=LSQCURVEFIT(FUN,X0,XDATA,YDATA) starts at X0 and finds coefficients X
%   to  best fit the nonlinear functions in FUN to the data YDATA (in the 
%   least-squares sense). FUN accepts inputs X and XDATA and returns a
%   vector (or matrix) of function values F, where F is the same size as
%   YDATA, evaluated at X and XDATA. NOTE: FUN should return FUN(X,XDATA)
%   and not the sum-of-squares sum((FUN(X,XDATA)-YDATA).^2).
%   ((FUN(X,XDATA)-YDATA) is squared and summed implicitly in the
%   algorithm.) 
%
%   X=LSQCURVEFIT(FUN,X0,XDATA,YDATA,LB,UB) defines a set of lower and
%   upper bounds on the design variables, X, so that the solution is in the
%   range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X=LSQCURVEFIT(FUN,X0,XDATA,YDATA,LB,UB,OPTIONS) minimizes with the
%   default parameters replaced by values in the structure OPTIONS, an
%   argument created with the OPTIMSET function. See OPTIMSET for details.
%   Used options are Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%   FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType,
%   LevenbergMarquardt, MaxFunEvals, MaxIter, DiffMinChange and
%   DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG,
%   PlotFcns, OutputFcn, and TypicalX. Use the Jacobian option to specify
%   that FUN also returns a second output argument J that is the Jacobian
%   matrix at the point X. If FUN returns a vector F of m components when X
%   has length n, then J is an m-by-n matrix where J(i,j) is the partial
%   derivative of F(i) with respect to x(j). (Note that the Jacobian J is
%   the transpose of the gradient of F.)
%
%   X = LSQCURVEFIT(PROBLEM) solves the non-linear least squares problem defined
%   in PROBLEM. PROBLEM is a structure with the function FUN in PROBLEM.objective,
%   the start point in PROBLEM.x0, the 'xdata' in PROBLEM.xdata, the 'ydata' in
%   PROBLEM.ydata, the lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub,
%   the options structure in PROBLEM.options, and solver name 'lsqcurvefit' in 
%   PROBLEM.solver. Use this syntax to solve at the command line a problem 
%   exported from OPTIMTOOL. The structure PROBLEM must have all the fields.
%
%   [X,RESNORM]=LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) returns the value of the
%   squared 2-norm of the residual at X: sum {(FUN(X,XDATA)-YDATA).^2}.
%
%   [X,RESNORM,RESIDUAL]=LSQCURVEFIT(FUN,X0,...) returns the value of residual,
%   FUN(X,XDATA)-YDATA, at the solution X. 
%
%   [X,RESNORM,RESIDUAL,EXITFLAG]=LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) returns 
%   an EXITFLAG that describes the exit condition of LSQCURVEFIT. Possible 
%   values of EXITFLAG and the corresponding exit conditions are
%
%     1  LSQCURVEFIT converged to a solution X.
%     2  Change in X smaller than the specified tolerance.
%     3  Change in the residual smaller than the specified tolerance.
%     4  Magnitude of search direction smaller than the specified tolerance.
%     0  Maximum number of function evaluations or of iterations reached.
%    -1  Algorithm terminated by the output function.
%    -2  Bounds are inconsistent.
%    -4  Line search cannot sufficiently decrease the residual along the
%         current search direction. 
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) 
%   returns a structure OUTPUT with the number of iterations taken in 
%   OUTPUT.iterations, the number of function evaluations in OUTPUT.funcCount, 
%   the algorithm used in OUTPUT.algorithm, the number of CG iterations (if 
%   used) in OUTPUT.cgiterations, the first-order optimality (if used) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA]=LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) 
%   returns the set of Lagrangian multipliers, LAMBDA, at the solution: 
%   LAMBDA.lower for LB and LAMBDA.upper for UB.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=LSQCURVEFIT(FUN,X0,XDATA,YDATA,...)
%   returns the Jacobian of FUN at X.
%
%   Examples
%     FUN can be specified using @:
%        xdata = [5;4;6];          % example xdata
%        ydata = 3*sin([5;4;6])+6; % example ydata
%        x = lsqcurvefit(@myfun, [2 7], xdata, ydata)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x,xdata)
%       F = x(1)*sin(xdata)+x(2);
%
%   FUN can also be an anonymous function:
%       x = lsqcurvefit(@(x,xdata) x(1)*sin(xdata)+x(2),[2 7],xdata,ydata)
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the curve-fitting 
%   problem given in the function myfun, which is parameterized by its second 
%   argument c. Here myfun is an M-file function such as
%
%       function F = myfun(x,xdata,c)
%       F = x(1)*exp(c*xdata)+x(2);
%
%   To solve the curve-fitting problem for a specific value of c, first assign 
%   the value to c. Then create a two-argument anonymous function that captures
%   that value of c and calls myfun with three arguments. Finally, pass this 
%   anonymous function to LSQCURVEFIT:
%
%       xdata = [3; 1; 4];           % example xdata
%       ydata = 6*exp(-1.5*xdata)+3; % example ydata     
%       c = -1.5;                    % define parameter 
%       x = lsqcurvefit(@(x,xdata) myfun(x,xdata,c),[5;1],xdata,ydata)
%
%   See also OPTIMSET, LSQNONLIN, FSOLVE, @, INLINE.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:27 $

%   The default algorithm when OPTIONS.LargeScale = 'off' is the 
%   Levenberg-Marquardt method with a mixed quadratic and cubic line search procedure.  
%   A Gauss-Newton method is selected by setting OPTIONS.LargeScale='off' and 
%   OPTIONS.LevenbergMarquardt='off'. 
%

% ------------Initialization----------------
defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off','FunValCheck','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','on', ...
   'OutputFcn',[],'PlotFcns',[]); 
% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   x = defaultopt;
   return
end

 if nargin < 7, options= [];
   if nargin < 6, UB = [];
      if nargin < 5, LB = [];
         end, end, end
 
problemInput = false;
 if nargin == 1
    if isa(FUN,'struct')
        problemInput = true;
        [FUN,x,XDATA,YDATA,LB,UB,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error('optim:lsqcurvefit:InputArg','The input to LSQCURVEFIT should be either a structure with valid fields or consist of at least four arguments.');
    end
end

if nargin < 4 && ~problemInput
    error('optim:lsqcurvefit:NotEnoughInputs', ...
        'LSQCURVEFIT requires four input arguments.');
end
% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(x,YDATA,LB,UB);
    if ~isequal('double', dataType)
        error('optim:lsqcurvefit:NonDoubleInput', ...
            'LSQCURVEFIT only accepts inputs of data type double.')
    end
catch
    error('optim:lsqcurvefit:NonDoubleInput', ...
        'LSQCURVEFIT only accepts inputs of data type double.')
end

if nargout > 5
   computeLambda = 1;
else 
   computeLambda = 0;
end

caller = 'lsqcurvefit';
[x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = ...
   lsqncommon(FUN,x,XDATA,YDATA,LB,UB,options,defaultopt,caller,...
              computeLambda,varargin{:});

