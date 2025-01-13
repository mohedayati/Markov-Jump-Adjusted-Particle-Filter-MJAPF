function state_jumps = jumping_strategy

%% Model details
dt=0.1;    tend = 100;    t=0:dt:tend;

%% Loading the pre-made transition probability matrix
load('C:\Users\hedayat2\Desktop\Simulation\Mohammad_Markovian\transition_probability_matrix.mat'); % Loads the transition probability matrix
%load('C:\Users\hedayat2\Desktop\Simulation\Mohammad_Markovian\transition_probability_matrix_train_CNN.mat'); % Loads the transition probability matrix
P = transition_prob;

%% Simulating the state jumps

%for i=1:1000

% Initialize the model with the beginning state (Nominal)
currentState = 1; % Assuming the system starts in state 1

N = length(t); % Number of time steps
states = zeros(1, N); % To record the state at each time step
states(1) = currentState;

for t = 2:N
    transitionProbabilities = P(currentState, :);
    nextState = randsample(1:9, 1, true, transitionProbabilities);
    states(t) = nextState;
    currentState = nextState;
end


% hold on
% legend
%end
state_jumps = states;