function mu = ky_spd_mean_chol(spd3d)

% See also KY_SPD_MEAN.


%% initialize
[p,~,N] = size(spd3d);
tmp_mat = zeros(p,p);

%% iterate
for n=1:N
    tmp_mat = tmp_mat + (1/N)*chol(spd3d(:,:,n));
end

%% finalize
mu = tmp_mat'*tmp_mat;

end