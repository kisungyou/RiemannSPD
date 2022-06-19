function logvecs = ky_spd_tangentialize(spd3d, Cref)

% KY_SPD_TANGENTIALIZE performs tangentialization of the data given a
% reference matrix. Note that the half-vectorization is applied for saving
% memory.
%
%   * USAGE
%       logvecs = KY_SPD_TANGENTIALIZE(spd3d)
%       logvecs = KY_SPD_TANGENTIALIZE(spd3d, Cref)
%
%   * INPUT
%       spd3d  a (p,p,N) 3d array of SPD matrices
%       Cref   a (p,p) SPD matrix. If not provided, identity is used.
%
%   * OUTPUT
%       logvecs a (N,p*(p+1)/2) matrix whose rows are tangential vectors.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [06/2022] initial implementation.


%% initialize
[p,~,N] = size(spd3d);
if (nargin < 2)
    Cref = eye(p);
end

%% iterate
logvecs = zeros(N,p*(p+1)/2);
for n=1:N
    tmp_vec = ky_spd_log(Cref, spd3d(:,:,n));
    logvecs(n,:) = ky_trf_vech(tmp_vec);
end


end