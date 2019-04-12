function [m,std] = krigeage(new_points,points,id)
% Compute mean and variance of a Gaussian process interpolating the
% funcion id = F(points)
%
% [m,std] = krigeage(points,id)
%
% Inputs: - new_points (m x d): points where to evaluate the process
%         - points (n x d): initial points where F has been computed
%         - id (n x 1) : values of F on points
%
%Outputs: - m (m x 1) mean of the Gaussian process on new_points
%         - std (m x 1) standard deviation on new_points
%

type='exp';
cL=3; %correlation length

M=covariance(points,points,cL,type);
N=covariance(new_points,points,cL,type);

m = N/M*id;
std = diag(1-N/M*N.');

end