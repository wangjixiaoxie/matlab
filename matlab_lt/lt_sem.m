%% LT 3/11/15 - find standard error of the mean for any-d array data
% 7/21/15 - MODIFIED to work with nans (will ignore, and downsample)
% WORKS FOR: 1d (any orientation), 2d (rows must be samples)

function Ysem=lt_sem(Y)
% Y is array or matrix. if matrix, then specify dimension to take using d;

% if ~exist('d','var');
%     d=1;
% end
Ysem=[];

if isempty(Y);
    Ysem=nan;
    
else
    if any(size(Y)==1);
        % is 1d
        % remove nans
        inds=isnan(Y);
        Y(inds)=[];
        
        % Continue;
        Ystd=std(Y);
        
        % get dimension
        sizes=size(Y);
        
        if sizes(1)==1;
            % then is wrong orientation
            Y=Y';
        end
        
        % get number of smaples (rows)
        N=size(Y);
        N=N(1);
        
        Ysem=Ystd./sqrt(N-1);
        
    else
        % then is 2d
        
        
        % do each columm separately
        numcols=size(Y,2);
        for i=1:numcols;
            Ycol=Y(:,i);
            
            % remove nans
            inds=isnan(Ycol);
            Ycol(inds)=[];
            
            % if nothing left, then were alll nans. put out a nan
            if isempty(Ycol);
                Ysem_tmp=nan;
            else
                Ystd=std(Ycol);
                
                N=length(Ycol);
                
                Ysem_tmp=Ystd./sqrt(N-1);
            end
            
            Ysem(i)=Ysem_tmp;
            
        end
    end
    
    
    
end
