function state_jumps = jumping_strategy

%% Model details
dt=0.1;    tend = 100;    t=0:dt:tend;

%% Loading the pre-made transition probability matrix
load('***.mat'); % Loads the transition probability matrix

P = transition_prob;

%% Simulating the state jumps

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

state_jumps = states;
