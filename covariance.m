function [ M ] = covariance( x, y, cor_length, type )
% Compute the covariance matrix
% [ M ] = covariance( x, y, cor_length, type )
%
% Inputs: - x : vector (n x d) of n points in dimension d
%         - y : vector (m x d) of m points in dimension d
%         - cor_length : correlation length
%         - type : kernel type (only 'exp' available)
%
% Output: - M (n x m) matrix of correlation

if ~strcmp(type,'exp')
    error('Kernel can only be ''exp'' !')
end

[n,~]=size(x);
[m,~]=size(y);

M=NaN(n,m);

for ix=1:n
   M(ix,:)=exp(-sqrt(sum((repmat(x(ix,:),m,1)-y) .^ 2,2))./cor_length);
end


end