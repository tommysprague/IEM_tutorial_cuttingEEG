function bb = build_basis_polar_mat(eval_at,chan_th,chan_size,chan_pow )
% adapted from build_basis_polar - this does everything at once, useful for
% generating a set of transformed/rotated basis functions

% will return bb, length(evalAt) x length(rfTh)
n_basis = length(chan_th);
[chan_th,eval_at] = meshgrid(squeeze(chan_th),squeeze(eval_at));

% TODO: make sure only one chan_size...


if nargin < 3 || isempty(chan_size)
    chan_size = 180; % default value
end

if nargin < 4 || isempty(chan_pow)
    chan_pow = n_basis-mod(n_basis,2);
end

% utility function to compute distance between two angles
ang_dist = @(a,b) min(mod(a-b,360),mod(b-a,360));


bb = (cosd( 180 * ang_dist(eval_at,chan_th) ./ (2*chan_size) ).^chan_pow) .* (ang_dist(eval_at,chan_th)<=chan_size) ;


return