function mu = ky_spd_mean_lerm(spd3d)

% See also KY_SPD_MEAN.


[p,~,N] = size(spd3d);
output = zeros(p,p);

for n=1:N
    output = output + real(logm(spd3d(:,:,n)))/N;
end

mu = real(expm(output));

end