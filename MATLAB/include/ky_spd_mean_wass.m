function mu = ky_spd_mean_wass(spd3d, maxiter, abstol)

% See also KY_SPD_MEAN.

% preproc
[p,~,N] = size(spd3d);
thr = double(abstol);
crit_inc  = 100000.00;
crit_iter = round(maxiter);


% run
mu_old = diag(diag(mean(spd3d,3)));

counter = 0;
while (crit_inc > thr)
    % compute the preliminary
    Skhalf = real(sqrtm(mu_old));
    Skhalfinv = inv(Skhalf);

    % summation part
    tmp_obj = zeros(p);
    for n=1:N
        tmp_obj = tmp_obj + real(sqrtm(Skhalf*spd3d(:,:,n)*Skhalf))/N;
    end
    
    % compute
    mu_new = Skhalfinv*tmp_obj*tmp_obj*Skhalfinv;

    % update
    crit_inc = norm(mu_old-mu_new,"fro");
    mu_old   = mu_new;
    counter = counter + 1;
    if (counter >= crit_iter)
        break;
    end
end

% return
mu = mu_old;

end