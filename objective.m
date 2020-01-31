function z_j = objective(x_j)
% The function "objective" return an objective vector z_j of x_j
% according  to a discrete stochastic nonlinear dynamical systems. 
% 
% z(t) = h(x(t)) for t = 0,...,T
% 
% Input: 
%   - x_j : state vector
% 
% Outputs : 
%   - z_j : objective vector 
% 
% Implemented example: 
%   z(t) = x(t) 
% 
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

z_j = x_j;
end