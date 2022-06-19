function dval = ky_wassdist(D, wx, wy, p)

% KY_WASSDIST computes the Wasserstein distance of order p given the cross
% pairwise distance matrix D.
%
%   * USAGE
%       DVAL = KY_WASSDIST(D)
%       DVAL = KY_WASSDIST(D, WX, WY)
%       DVAL = KY_WASSDIST(D, WX, WY, P)
%
%   * INPUT
%       D      an (m,n) cross pairwise distance.
%       wx     (optional) a length-m vector of weights. Default is uniform.
%       wy     (optional) a length-n vector of weights. Default is uniform.
%       p      (optional) the order of Wasserstein distance.
%
%   * OUTPUT
%       dval   the p-Wasserstein distance.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [06/2022] initial implementation.


%% initialize
dxy = D; % pdist2, probably
[m,n] = size(D);

if (nargin < 2)
    wx = ones(m,1)/m;
end
if (nargin < 3)
    wy = ones(n,1)/n;
end

wx  = wx/sum(wx); wx = wx(:); % vector is a single-column matrix
wy  = wy/sum(wy); wy = wy(:); 

if (nargin < 4)
    p = 2.0;
end

cxy   = dxy.^p;

c  = cxy(:);
A1 = kron(ones(1,n), eye(m));
A2 = kron(eye(n), ones(1,m));
A  = [A1;A2];

f_obj = c;
f_con = A;
f_rhs = [wx;wy]; f_rhs = f_rhs(:);
f_lwb = zeros(length(f_obj),1);
f_upb = ones(length(f_obj),1);
f_opt = optimoptions("linprog","Display","none");
f_sol = linprog(f_obj, [], [], f_con, f_rhs, f_lwb, f_upb, f_opt);

dval  = (sum(f_sol.*c))^(1/p);


end