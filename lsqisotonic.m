function yhat = lsqisotonic(x,y,w)
%LSQISOTONIC Isotonic least squares.
%   YHAT = LSQISOTONIC(X,Y) returns a vector of values that minimize the
%   sum of squares (Y - YHAT).^2 under the monotonicity constraint that
%   X(I) > X(J) => YHAT(I) >= YHAT(J), i.e., the values in YHAT are
%   monotonically non-decreasing with respect to X (sometimes referred
%   to as "weak monotonicity").  LSQISOTONIC uses the "pool adjacent
%   violators" algorithm.
%
%   If X(I) == X(J), then YHAT(I) may be <, ==, or > YHAT(J) (sometimes
%   referred to as the "primary approach").  If ties do occur in X, a plot
%   of YHAT vs. X may appear to be non-monotonic at those points.  In fact,
%   the above monotonicity constraint is not violated, and a reordering
%   within each group of ties, by ascending YHAT, will produce the desired
%   appearance in the plot.
%
%   YHAT = LSQISOTONIC(X,Y,W) performs weighted isotonic regression using
%   the non-negative weights in W.

%   Copyright 2003-2006 The MathWorks, Inc.


%   References:
%      [1] Kruskal, J.B. (1964) "Nonmetric multidimensional scaling: a
%          numerical method", Psychometrika 29:115-129.
%      [2] Cox, R.F. and Cox, M.A.A. (1994) Multidimensional Scaling,
%          Chapman&Hall.

n = numel(x);
if nargin<3
    yclass = superiorfloat(x,y);    % 确定数据格式是double还是single
else
    yclass = superiorfloat(x,y,w);
end

% Sort points ascending in x, break ties with y.
% 将x的点升序排列，打断与y的联系 
[xyord,ord] = sortrows([x(:) y(:)]); % 将x,y拼在一起进行排列，每列都升序排列
iord(ord) = 1:n;    % 记录排列的顺序
xyord = double(xyord);

% Initialize fitted values to the given values.
% 用给定的值来初始化匹配的值，只用前两列？
yhat = xyord(:,2);

block = 1:n;
if (nargin == 3) && ~isempty(w)
    w = double(w(:)); w = w(ord); % reorder w as a column

    % Merge zero-weight points with preceding pos-weighted point (or
    % with the following pos-weighted point if at start).
    posWgts = (w > 0);
    if any(~posWgts)
        idx = cumsum(posWgts); idx(idx == 0) = 1;
        w = w(posWgts);
        yhat = yhat(posWgts);
        block = idx(block);
    end

else
    w = ones(size(yhat));
end

while true
    % If all blocks are monotonic, then we're done.
    diffs = diff(yhat); %求行间差
    if all(diffs >= 0), break; end

    % Otherwise, merge blocks of non-increasing fitted values, and set the
    % fitted value within each block equal to a constant, the weighted mean
    % of values in that block.
    idx = cumsum([1; (diffs>0)]);
    sumyhat = accumarray(idx,w.*yhat);
    w = accumarray(idx,w);
    yhat = sumyhat ./ w;
    block = idx(block);
end

% Broadcast merged blocks out to original points, and put back in the
% original order and shape.
% 将融合的数据排回原来的位置
yhat = yhat(block);
yhat = reshape(yhat(iord), size(y));
if isequal(yclass,'single')
    yhat = single(yhat);
end
