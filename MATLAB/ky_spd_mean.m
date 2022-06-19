function mu = ky_spd_mean(spd3d, varargin)

% KY_SPD_MEAN is a unifying function to compute the mean of a collection of
% SPD-valued matrices under a number of geometries that are listed below.
% 
% computes the geometric mean of SPD matrices under a number of
% geometries
% 2-Wasserstein metric, which is a pseudo-Riemannian geometry. The
% algorithm implemented here is the one proposed by 
%
%   * USAGE
%       mu = KY_SPD_MEAN(spd3d)
%       mu = KY_SPD_MEAN(spd3d, 'PARAM1',val1, 'PARAM2',val2, ...)
%
%   * INPUT
%       spd3d  a (p,p,N) 3d array of SPD matrices.
%
%   * OPTIONS (Parameter-Value pairs)
%       'geometry' - the name of corresponding geometry. We allow
%         case-insensitive input from one of the following choices:
%           'airm'        - affine-invariant Riemannian metric.
%           'chol'        - Cholesky decomposition.
%           'euclidean'   - (default) Euclidean geometry.
%           'harmonic'    - harmonic mean.
%           'lerm'        - log-Euclidean Riemannian metric.
%           'wasserstein' - 2-Wasserstein geometry. The algorithm
%                           implemented is by Alvarez-Esteban et al (2016).
%
%       'maxiter' - the maximum number of iterations. (Default: 100).
%
%       'abstol' - the stopping criterion. (Default: 1e-8).
%
%   * OUTPUT
%       mu     a (p,p) mean matrix under the geometry.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [06/2022] initial implementation.



%% initialize
%   default options
par_geometry = char("euclidean");
par_maxiter  = 100;
par_abstol   = 1e-8;

%   varargin
for setting=1:2:length(varargin)
    switch varargin{setting}
        case 'geometry'
            par_geometry = char(lower(varargin{setting + 1}));
        case 'maxiter'
            par_maxiter = round(varargin{setting+1});
        case 'abstol'
            par_abstol = double(varargin{setting+1});
    end
end

%   adjust absurd inputs
if (par_maxiter < 2)
    par_maxiter = 5;
end
if ((par_abstol<=0)||(par_abstol>1))
    par_abstol = 1e-8;
end

%% run per geometry
switch par_geometry
    case "airm"
        mu = ky_spd_mean_airm(spd3d, par_maxiter, par_abstol);
    case "chol"
        mu = ky_spd_mean_chol(spd3d);
    case "euclidean"
        mu = ky_spd_mean_euclid(spd3d);
    case "harmonic"
        mu = ky_spd_mean_harmonic(spd3d);
    case "lerm"
        mu = ky_spd_mean_lerm(spd3d);
    case "wasserstein"
        mu = ky_spd_mean_wass(spd3d, par_maxiter, par_abstol);
    otherwise
        error("* ky_spd_mean : no such geometry is available.")
end
end