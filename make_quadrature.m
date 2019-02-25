function [ points, weights ] = make_quadrature( dim, order, tensorisation )
% [ points, weights ] = make_quadrature( dim, order(,tensorisation) )
%
%   Use the Tasmanian package to generate a Gauss Hermite quadrature
%
%   dim: dimension
%   order:  max order retained
%   tensorisation can be: -'tensor'=full (default)
%                         -'qptotal','level','curved','tensor','iptensor',
%                          'iptotal','ipcurved','qptotal','qpcurved',
%                          'hyperbolic','iphyperbolic','qphyperbolic'=sparse
% 

if nargin<3
    tensorisation='tensor';
end

prec = order;
domain = repmat([ 0, 1/2],[dim,1]); %used to have exp(-x^2/2) for prob and exp(-x^2) for phys
ordr = 0; % not used by Gauss-hermite rule
[ weights, points ] = tsgMakeQuadrature( dim, 'gauss-hermite',tensorisation, prec, ordr,domain,0);
weights=weights./((2*pi)^(dim/2)); %normalization pi for phys 2*pi for prob

end

