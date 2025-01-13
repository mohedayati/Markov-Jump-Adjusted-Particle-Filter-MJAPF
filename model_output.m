function [newState, newState_w] = model_output(t_inquiry, mode_at_t, particles_history, particles_history_w, modeHistory, p)

%states = jumping_strategy;
states_nominal = ones(1,1001);
states = ones(1,1001);
mode_p_till_t_minus_1 = modeHistory(1:t_inquiry-1, p);
states(1:t_inquiry-1) = mode_p_till_t_minus_1;
states(t_inquiry) = mode_at_t;



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

%% Modeling the process noise


A_process_1 = 0;
A_measurement_1 = 0;

B_process_1 = 1.5*0.05;
B_measurement_1 = 30*0.05;

name = 'Logistic'; % Type of distribution


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


ti = (t_inquiry-1)*dt-2*dt;
u = 5*sin(0.2*ti);
x(:,i+1) = f_model(x(:,i),wf(:,i),u,dt); % + er_process; % wf is responsbile for injecting the faults | x is the combined state vector for I_dot and omega_dot
% ************* process error has been omitted from the output of the model_output.m function


X_simulated = x(1,:); % only I
X_simulated_w = x(2,:); % only w
newState = X_simulated(t_inquiry);
newState_w = X_simulated_w(t_inquiry);
