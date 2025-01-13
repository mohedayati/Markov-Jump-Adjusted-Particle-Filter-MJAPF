function [newState, newState_w] = model_output(t_inquiry, mode_at_t, particles_history, particles_history_w, modeHistory, p)

%states = jumping_strategy;
states_nominal = ones(1,1001);
states = ones(1,1001);
mode_p_till_t_minus_1 = modeHistory(1:t_inquiry-1, p);
states(1:t_inquiry-1) = mode_p_till_t_minus_1;
states(t_inquiry) = mode_at_t;

% for h = 1:length(states)
%     if states(h) ~= 1
%         break
%     end
% end


n = 2;                                         %number of state
m = 2;                                         %numer of measurements
w(:,1)     = [.01;2];
dt=0.1;    tend = 100;    t=0:dt:tend;
tsim = 0:dt:tend-dt;
phi=10;
r = 1e-4;                                      %std of measurement

er = zeros(n,length(tsim));
er_process = zeros(n,length(tsim));
x = zeros(n,length(tsim)+1);
uc = zeros(3,length(tsim));
z = zeros(n,length(tsim));
wf = ones(length(w),length(tsim));
wf(:,1)=[0.029, 8];

% %% Making a fully nominal model for optimization in generating faulty outputs
% 
% for i=1:length(states_nominal)
%     if i~=1 & states_nominal(i) == states_nominal(i-1)
%         continue
%     else
%         [ftime, wd]= faults(i/10,tend+dt, states_nominal(i));
%         for kk=1:length(ftime)-1
%             k=find(t>=ftime(kk) & t<=ftime(kk+1));
%             for j=k
%                 wf(:,j)=wd(:,kk);
%             end
%         end
%     end
% end
% 
% clear wf wd ftime
% 
% wf = ones(length(w),length(tsim));
% wf(:,1)=[0.029, 8];

%% Modeling the process noise

% Laplace (Double Exponential) noise:
% Laplace noise parameters for two variables
% Process noise parameters
mu_process1 = 0; sigma_process1 = 50*r*phi*1e-1*10; % Mean and SD for I_dot
mu_process2 = 0; sigma_process2 = 50*r*phi*1e-1*10; % Mean and SD for omega_dot

A_process_1 = 0;
A_measurement_1 = 0;

B_process_1 = 1.5*0.05;
B_measurement_1 = 30*0.05;

name = 'Logistic'; % Type of distribution



%er_process = random(name, A_process_1, B_process_1); % Additive Logistic process noise for both variables


%% Generating faulty model outputs for each particle


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

tsim_truncated = tsim(1:t_inquiry-1);


x(:,1)=[0;0]; 
x(1,1:t_inquiry-1)= particles_history(1:t_inquiry-1, p);
x(2,1:t_inquiry-1)= particles_history_w(1:t_inquiry-1, p);

i = t_inquiry-1;


%    i=i+1;
    % Gaussian noise:
    %er_process = 50*[normrnd(0,r*phi*1e-1*10); normrnd(0,r*phi*1e-1*10)]; % process noise

    ti = (t_inquiry-1)*dt-2*dt;
    u = 5*sin(0.2*ti);
    %uc(:,i+1) = u;
    %wf
    x(:,i+1) = f_model(x(:,i),wf(:,i),u,dt); % + er_process; % wf is responsbile for injecting the faults | x is the combined state vector for I_dot and omega_dot
    % ************* process error has been omitted from the output of the
    % model_output.m function
    %z(:,i+1) = x(:,i+1) + er(:,i+1);


% for ti=tsim_truncated-2*dt
%     i=i+1;
%     % Gaussian noise:
%     er_process(:,i+1) = 50*[normrnd(0,r*phi*1e-1*10); normrnd(0,r*phi*1e-1*10)]; % process noise
%     er(:,i+1) = 50*[normrnd(0,r*phi*1e-1); normrnd(0,r*phi*1e-1)]; % measurement noise
% 
%     u = 5*sin(0.2*ti);
%     uc(:,i+1) = u;
%     x(:,i+1) = f_model(x(:,i),wf(:,i),u,dt) + er_process(:,i+1); % wf is responsbile for injecting the faults | x is the combined state vector for I_dot and omega_dot
%     %z(:,i+1) = x(:,i+1) + er(:,i+1);
% end

% X=x;
X_simulated = x(1,:); % only I
X_simulated_w = x(2,:); % only w
newState = X_simulated(t_inquiry);
newState_w = X_simulated_w(t_inquiry);
clear X_simulated states
% Z=z;
% Z_simulated=z(1,:); % only I_dot
