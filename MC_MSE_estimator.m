function [mse] = MC_MSE_estimator(meas,T,n_draw,n_part)
% Monte Carlo estimator of the expected mean square error of a particle filter 
% estimate of the state vector of a stochastic nonlinear dynamical systems
% with a given set of measurements.
%
% A discrete stochastic nonlinear dynamical system is modelled:
%   x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
%   y(t) = g(x(t),v(t)) for t in a given measurment set
%   z(t) = h(x(t)) for t = 0,...,T 
%   x(0) ~ F 
%   w(t) and v(t) are the process and measurement noise respectively witn
%   known probability density functions
% 
% The stochastic nonlinear dynamical systems is specified by the functions: 
%   - model(T) : draw a state vector x_j 
%   - measurements(x_j,T) : draw a measurement vector y_j of x_j
%   - objective(x_j) : return an objective z_j of x_j
%
% n_draw pairs of random variables (x_j,y_j) are simulated according to the
% stochastic nonlinear dynamical model. 
% For each draw (x_j,y_j), an estimate of x_j, tau_j, is computed using a 
% particle filter with N particles. 
% 
% The particle filter is specified by the function: 
%   - particle_filter(y_j,measurement_times,n_part,T)
%
% The MSE of the particle filter is for each draw. 
% The estimate of the expected MSE is computed by averaging the mse over
% the n_draw draws. 
%
% Inputs: 
%   - meas : index of measurements
%   - T : length of the time interval 
%   - n_draw : number of draws
%   - n_part : number of particles used in the particle filter
% 
% Output: 
%   - mse : estimation of the expected mean square error
%
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

measurement_times = zeros(1,T+1); 
measurement_times(meas+1) = 1;

mse = 0;
for j = 1:n_draw
        %1.Simulation

        %1.1. Random motion model : X_j
        x_j = model(T);

        %1.2. Artificial data record Y_j
        y_j = measurements(x_j,T);

        %2. Filtering
        tau_j = particle_filter(y_j,measurement_times,n_part,T);

        %3. MSE computation
        if isnan(tau_j)
            mse = nan;%10^8;
            disp('All particles have w = 0')
            break
        else
            mse = mse + 1/n_draw*mean((objective(x_j)-objective(tau_j)).^2,2);          
        end
end
end
