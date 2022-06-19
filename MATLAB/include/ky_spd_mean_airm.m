function mu = ky_spd_mean_airm(spd3d, maxiter, abstol)

% See also KY_SPD_MEAN.


%% initialize
stop_thr  = double(abstol);
stop_inc  = 10000.0;
stop_iter = round(maxiter);
[p,~,N] = size(spd3d);

mu_old = diag(diag(mean(spd3d, 3)));

%% iterate
counter = 0;
while (stop_inc > stop_thr)
    mu_tmp = zeros(p);
    for n=1:N
        mu_tmp = mu_tmp + (1/N)*ky_spd_log(mu_old, spd3d(:,:,n));
    end
    mu_new = ky_spd_exp(mu_old, mu_tmp);

    stop_inc = norm(mu_new - mu_old,"fro");
    mu_old   = mu_new;

    counter = counter + 1;
    if (counter >= stop_iter)
        break;
    end
end

mu = mu_old;
end
