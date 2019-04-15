function [um, um_std, gPC_mse] = gPC_mod_LS(u,Psi,varargin)
%GPC_MOD_LS returns a set of spectral gPC modal coefficients from a given realization ensemble
%associated to a continuous probability measure via a Least-Square regression.
%It also returns estimates of the std errors for those coefficients,
%and an estimate of the std of the regression error term
%   GPC_MOD_LS(U,PSI,VARARGIN), returns values of U in the modal space: Um,
%   std errors: um_std and mean square error: gPC_mse
%   The first argument are the realizations of the RV(s) U. Each row of U is a
%   realization and each column a RV.
%   The second argument represents the gPC polynomial basis evaluated at some points
%   (realizations of the underlying RV(s) seed(s)). Each row of Psi is a realization and
%   each column a polynomial member.
%
%   Remark: vector/matrice arguments must be properly ordered
weight_flag = false;
if nargin == 3,
    weight = varargin{:}; weight_flag = true;
end

[Nq,P] = size(Psi); [Nq,n]  = size(u);

um = zeros(P,n); um_std = um; gPC_mse = zeros(1,n);

% um = Psi\u;
if weight_flag
    for i=1:n,
        [um(:,i),um_std(:,i), gPC_mse(i)] = lscov(Psi,u(:,i),weight);
    end
else
    [um,um_std, gPC_mse] = lscov(Psi,u);
end
