clear
clc

%%%%% READ ME IMPORTANTTTT

% I HAVE MODELED THE MARKOVIAN JUMPS AS PROCESS NOISE TOO. THE REASON BEING
% THAT MARKOVIAN JUMPS *DON'T PROPAGATE THROUGH TIME* DUE TO ONLY DEPENDING ON
% THE PREVIOUS TIME-STEP. WE DON'T INCLUDE THE PROPAGATED
% PROCESS NOISE IN THE PRIOR DISTRIBUTION ( p(x_i given x_(i-1)) ) BUT WE
% INCLUDE THE ADDED PROCESS NOISE AT EACH TIME STEP IN THE PRIOR. WE TREAT
% THE MARKOVIAN JUMP TRANSITIONS LIKE THE ADDED PROCESS NOISE.
%% No process noise

n = 2;                                         %number of state
m = 2;                                         %numer of measurements
w(:,1)     = [.01;2];
dt=0.1;    tend = 100;    t=0:dt:tend;
tsim = 0:dt:tend-dt;
T = length(t);
phi=10;
r = 1e-4;                                      %std of measurement

er = zeros(n,length(tsim));
er_process = zeros(n,length(tsim));
x = zeros(n,length(tsim));
uc = zeros(3,length(tsim));
z = zeros(n,length(tsim));
wf = ones(length(w),length(tsim));
wf(:,1)=[0.029, 8];



[ftime, wd]= faults(0,tend, 1);
for kk=1:length(ftime)-1
    k=find(t>=ftime(kk) & t<=ftime(kk+1));
    for j=k
        wf(:,j)=wd(:,kk);
    end
end


i=0; x(:,1)=[0;0];
for ti=tsim-2*dt
    i=i+1;
    u = 5*sin(0.2*ti);
    uc(:,i+1) = u;
    x(:,i+1) = f_model(x(:,i),wf(:,i),u,dt); % wf is responsbile for injecting the faults | x is the combined state vector for I_dot and omega_dot
end

% X=x;
X_noise_less=x(1,:); % only I
X_noise_less_w=x(2,:); % only w
% figure
% plot(X_noise_less)

%% Simulation of a model for calculating the KDE

clear x z uc er er_process wf

num_seq = 10000;

residuals_all = zeros(num_seq, 1001);

    n = 2;                                         %number of state
    m = 2;                                         %numer of measurements
    w(:,1)     = [.01;2];
    dt=0.1;    tend = 100;    t=0:dt:tend;
    tsim = 0:dt:tend-dt;
    T = length(t);
    phi=10;
    r = 1e-4;                                      %std of measurement

for h=1:num_seq
    states = jumping_strategy_v;


    er = zeros(n,length(tsim));
    er_process = zeros(n,length(tsim));
    x = zeros(n,length(tsim));
    uc = zeros(3,length(tsim));
    z = zeros(n,length(tsim));
    wf = ones(length(w),length(tsim));
    wf(:,1)=[0.029, 8];



    for i=1:length(states)
        if i~=1 & states(i) == states(i-1)
            continue
        else
            [ftime, wd]= faults(i/10,tend+dt, states(i));
            for kk=1:length(ftime)-1
                k=find(t>=ftime(kk) & t<=ftime(kk+1));
                for j=k
                    wf(:,j)=wd(:,kk);
                end
            end
        end
    end


    % Process noise parameters
    mu_process1 = 0; sigma_process1 = 50*r*phi*1e-1*10; % Mean and SD for I_dot
    mu_process2 = 0; sigma_process2 = 50*r*phi*1e-1*10; % Mean and SD for omega_dot

    A_process_1 = 0;

    B_process_1 = 2*0.1*0.05;






    name = 'Logistic'; % Type of distribution


    i=0; x(:,1)=[0;0];
    for ti=tsim-2*dt
        i=i+1;
        er_process(:,i+1) = [random(name, A_process_1, B_process_1); random(name, A_process_1, 10*B_process_1)]; % Additive Laplace process noise for both variables
        u = 5*sin(0.2*ti);
        uc(:,i+1) = u;
        x(:,i+1) = f_model(x(:,i),wf(:,i),u,dt) + er_process(:,i+1); % wf is responsbile for injecting the faults | x is the combined state vector for I_dot and omega_dot

    end

    % X=x;
    X_noisy=x(1,:); % only I
    X_noisy_w=x(2,:); % only w

    residual = X_noisy - X_noise_less;
    residuals_all(h,:) = residual;

    residual_w = X_noisy_w - X_noise_less_w;
    residuals_all_w(h,:) = residual_w;
    h
end

% hold on
% plot(X_noisy)


data = reshape(residuals_all, 1, []);
data_w = reshape(residuals_all_w, 1, []);


% Histogram for empirical PDF (I)
figure;
histogram(data, 'Normalization', 'pdf', 'BinWidth', 0.001); % 'pdf' normalizes the histogram
hold on;

% Kernel Density Estimation for a smoother PDF
% Define the specific points 'xi' where we want to evaluate the KDE
xi = linspace(min(data), max(data), 100000); % 1000 evenly spaced points over the range
[f, xi] = ksdensity(data);
plot(xi, f, 'LineWidth', 2);

title('Empirical PDF (I)');
xlabel('Measurements');
ylabel('Probability Density');
legend('Histogram', 'KDE');
hold off;

% Histogram for empirical PDF (w)
figure;
histogram(data_w, 'Normalization', 'pdf', 'BinWidth', 0.001); % 'pdf' normalizes the histogram
hold on;

% Kernel Density Estimation for a smoother PDF
% Define the specific points 'xi' where we want to evaluate the KDE
xi_w = linspace(min(data_w), max(data_w), 100000); % 1000 evenly spaced points over the range
[f_w, xi_w] = ksdensity(data_w);
plot(xi_w, f_w, 'LineWidth', 2);

title('Empirical PDF (w)');
xlabel('Measurements');
ylabel('Probability Density');
legend('Histogram', 'KDE');
hold off;

integral_kde = trapz(xi, f)
integral_kde_w = trapz(xi_w, f_w)

kde_pdf_function = @(x) interp1(xi, f, x, 'linear', 'extrap');
kde_pdf_function_w = @(x) interp1(xi_w, f_w, x, 'linear', 'extrap');

x_range = linspace(min(data) - 1, max(data) + 1, 1000);  % Extend the range slightly beyond the min and max of the data
kde_values = kde_pdf_function(x_range);

x_range_w = linspace(min(data_w) - 1, max(data_w) + 1, 1000);  % Extend the range slightly beyond the min and max of the data
kde_values_w = kde_pdf_function_w(x_range_w);

figure;  % Create a new figure window
plot(x_range, kde_values, 'LineWidth', 2);  % Plot the KDE
title('Kernel Density Estimate (I)');
xlabel('Data Values');
ylabel('Density');
grid on;  % Turn on the grid


figure;  % Create a new figure window
plot(x_range_w, kde_values_w, 'LineWidth', 2);  % Plot the KDE
title('Kernel Density Estimate (w)');
xlabel('Data Values');
ylabel('Density');
grid on;  % Turn on the grid

% Save the function and necessary variables to a MAT-file
save('kdeFunction.mat', 'kde_pdf_function');
save('kdeFunction_w.mat', 'kde_pdf_function_w');