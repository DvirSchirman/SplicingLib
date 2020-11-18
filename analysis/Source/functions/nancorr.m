function [r,p] = nancorr(x,y,type)

if nargin<3
    type='Pearson';
end

if nargin==1
    [r,p]=corr(x,'rows','pairwise','type',type);
elseif nargin==2
    if isstr(y)
        type=y;
        [r,p]=corr(x,'rows','pairwise','type',type);
    else
        [r,p]=corr(x,y,'rows','pairwise','type',type);
    end
elseif nargin==3
    [r,p]=corr(x,y,'rows','pairwise','type',type);
end