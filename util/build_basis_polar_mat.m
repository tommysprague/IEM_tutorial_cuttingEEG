function bb = build_basis_polar_mat(evalAt,rfTh)
% adapted from build_basis_polar - this does everything at once, useful for
% generating a set of transformed/rotated basis functions

% will return bb, length(evalAt) x length(rfTh)
n_basis = length(rfTh);
[rfTh,evalAt] = meshgrid(squeeze(rfTh),squeeze(evalAt));

%bb = nan(length(rfTh),1);
%for ii = 1:length(rfTh)
    bb = cosd(0.5*min(mod(evalAt-rfTh,360),mod(rfTh-evalAt,360))).^(n_basis-mod(n_basis,2));
%end
return