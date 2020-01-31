function part = initialization(n_part)
% Function initializing the particles according to the model 
% X_0 ~ Xsi(x_0) dx_0
% 
% Input : 
%   - n_part : number of particles
% 
% Output : 
%   - part : particles
% 
% Date : 23/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

part = randn(1,n_part)*5 ; 

end