function mu = ky_spd_mean_harmonic(spd3d)

% See also KY_SPD_MEAN.


[p,~,N] = size(spd3d);
mu_inv = zeros(p);
for n=1:N
    mu_inv = mu_inv + inv(spd3d(:,:,n))/N;
end
mu = inv(mu_inv); % at the last step, don't forget!


end