%% Plot function 
close all;
clc;

n_draw_comp = 10^4; %number of draws to build the histogram
n_measurements = 21; %number of observations 
T = 60; %number of time steps 

n_part = 250; %number of particles
n_draw = 100; %number of draws for the MC MSE estimator 

pop_size = 50; %population size of the GA algorithm 
max_gen = 25;
n_eval = pop_size*max_gen;

%% 1. Computation of the measurement times reg, GA, RT 
display('1. Computation of the measurement times reg, GA, RT');

meas_reg = round(linspace(0,T,n_measurements));

display('GA computations started');
t_start_GA = tic;
[meas_GA,~,avgCostHist_GA ,minCostHist_GA] = genetical_algo(n_measurements,T,pop_size,max_gen,n_part,n_draw);
t_elapsed_GA = toc(t_start_GA);
display(['GA computations completed, time elapsed = ',num2str(t_elapsed_GA,'%.0f sec')]);

display('RT computations started');
t_start_RT = tic;
[meas_RT,~,avgCostHist_RT,minCostHist_RT] = random_trials(T,n_measurements, n_eval, n_part, n_draw);
t_elapsed_RT = toc(t_start_RT);
display(['RT computations completed, time elapsed = ',num2str(t_elapsed_RT,'%.0f sec')]);

%% 2. Comparison of the performances of the methods by running n_draw_comp draws
display(['2. Comparison of the performances of the methods by running ',num2str(n_draw_comp,'%.0f draws')]);

measurements_reg = zeros(1,T+1); 
measurements_reg(meas_reg+1) = 1;

measurements_GA = zeros(1,T+1); 
measurements_GA(meas_GA+1) = 1;

measurements_RT = zeros(1,T+1); 
measurements_RT(meas_RT+1) = 1;

%0. Constants definition 

mse_reg = 0;
mse_GA = 0; 
mse_RT = 0;

gains_GA = zeros(n_draw_comp,1);
gains_RT = zeros(n_draw_comp,1);

for j = 1:n_draw_comp
        %1.Simulation

        %1.1. Random motion model : X_j
        x_j = model(T);

        %1.2. Artificial data record Y_j
        y_j = measurements(x_j,T);

        %2. Filtering
        tau_j_reg = particle_filter(y_j,measurements_reg,n_part,T);
        tau_j_GA = particle_filter(y_j,measurements_GA,n_part,T);
        tau_j_RT = particle_filter(y_j,measurements_RT,n_part,T);
        
        %3. MSE computation
        err_reg = mean((objective(x_j)-objective(tau_j_reg)).^2,2);
        err_GA =  mean((objective(x_j)-objective(tau_j_GA)).^2,2);
        err_RT = mean((objective(x_j)-objective(tau_j_RT)).^2,2);
        
        mse_reg = mse_reg + 1/n_draw_comp*err_reg;
        mse_GA = mse_GA + 1/n_draw_comp*err_GA; 
        mse_RT = mse_RT + 1/n_draw_comp*err_RT;
        
        gains_GA(j) = (err_reg - err_GA)/err_reg;
        gains_RT(j) = (err_reg - err_RT)/err_reg;
end
  
%% 3. Plots
fontsize = 7*2;

figure;    
histogram(gains_GA,'Normalization','pdf')
positive_gain_GA = zeros(n_draw_comp,1);
positive_gain_GA(gains_GA>=0) = 1/n_draw_comp;
positive_gain_GA = sum(positive_gain_GA);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with GA algorithm')

display(['GA : average gain = ' num2str(mean(gains_GA)*100,'%.1f %%')]);
display(['GA : fraction of positive gain = ' num2str(positive_gain_GA*100,'%.1f %%')]);

figure;
histogram(gains_RT,'Normalization','pdf')
positive_gain_RT = zeros(n_draw_comp,1);
positive_gain_RT(gains_RT>=0) = 1/n_draw_comp;
positive_gain_RT = sum(positive_gain_RT);
line([0 0], get(gca, 'ylim'),'Color','red');
ylabel('probability density function','interpreter','latex')
xlabel('relative gain, $g$','interpreter','latex')
title('Histogram of the relative gain obtained with RT algorithm')

display(['RT : average gain = ' num2str(mean(gains_RT)*100,'%.1f %%')]);
display(['RT : fraction of positive gain = ' num2str(positive_gain_RT*100,'%.1f %%')]);


FIG = figure; hold on;
genList = (1:length(avgCostHist_GA))-1;
evalList = (1:n_eval);
plot(genList*pop_size,0*avgCostHist_GA+mse_reg,'-.k');
plot(genList*pop_size,avgCostHist_GA,'--r','LineWidth',1.5);
plot(genList*pop_size,minCostHist_GA,'-r','LineWidth',1.5);
plot(evalList,avgCostHist_RT,'--b');
plot(evalList,minCostHist_RT,'-b');
box on
xlim([0 n_eval])
legend('RT: average cost','RT: min cost','GA: average cost', 'GA: min cost GA', 'regualar MSE cost')
xlabel('number of cost function evaluations','interpreter','latex')
ylabel('$\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
title('Evolution of the average and minimum cost $\hat{\mathrm{E}}_{\mathrm{MSE}}$','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'Fontsize',fontsize);