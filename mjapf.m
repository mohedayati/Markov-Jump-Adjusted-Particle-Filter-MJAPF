function Z_filtered = mjapf(Z_real, Z_real_w, numParticles, kde_f_I, kde_f_W)

poolObj = gcp('nocreate'); % Gets current pool without creating a new one if it doesn't exist

cmap = colororder();
kde_pdf_function = kde_f_I;
kde_pdf_function_w = kde_f_W;
%% Simulation of a model for generating hypothetical measurements


states = jumping_strategy;

n = 2;                                         %number of state
m = 2;                                         %numer of measurements
w(:,1)     = [.01;2];
dt=0.1;    tend = 100;    t=0:dt:tend;
tsim = 0:dt:tend-dt;
T = length(t);
phi=10;
r = 1e-4;                                      %std of measurement



A_process_1 = 0;
A_measurement_1 = 0;

B_process_1 = 2*0.1*0.05;
B_measurement_1 = 4*0.05;

name = 'Logistic'; % Type of distribution



%% Forming the importance density function



data = zeros(1, 1000);
for xii = 1:1000
    data(xii) = (1-0.5)*random(name, A_process_1, 1*B_process_1) - 0.5*random(name, A_measurement_1, 1*B_measurement_1);
end

xi_kde = linspace(min(data), max(data), 1000); % 1000 evenly spaced points over the range
[f_kde, xi_kde] = ksdensity(data);

integral_kde_q = trapz(xi_kde, f_kde);

kde_pdf_function_q = @(x) interp1(xi_kde, f_kde, x, 'linear', 'extrap');


data_w = zeros(1, 1000);
for xii = 1:1000
    data_w(xii) = (1-0.5)*random(name, A_process_1, 1*10*B_process_1) - 0.5*random(name, A_measurement_1, 1*10*B_measurement_1);
end

xi_kde_w = linspace(min(data_w), max(data_w), 1000); % 1000 evenly spaced points over the range
[f_kde_w, xi_kde_w] = ksdensity(data_w);

integral_kde_q_w = trapz(xi_kde_w, f_kde_w);

kde_pdf_function_q_w = @(x) interp1(xi_kde_w, f_kde_w, x, 'linear', 'extrap');


%% Simple Jump Markov Particle Filter



% Define system parameters
% numParticles = 50;  % Number of particles
N_eff_ratio = 1;
N_eff_thresh = floor(numParticles * N_eff_ratio);
numModes = 9;  % Number of modes
load('transition_probability_matrix.mat'); % Loads the transition probability matrix
transitionMatrix = transition_prob;
measurementNoiseMean = 0;  % Mean of the measurement noise
measurementNoiseSD = 50*r*phi*1e-1;  % Standard deviation of the measurement noise

% Initialize particles
particles = rand(1, numParticles);  % Random initial state
particles_w = rand(1, numParticles);  % Random initial state
mode = randi([1, numModes], 1, numParticles);  % Initial mode for each particle


% Initialization with specific state range
initialStateRange = [-1*0.01, 1*0.01];  % Initial state range (I)
initialStateRange_w = [-6*0.01, 6*0.01];  % Initial state range (omega or w)

% Initialize particles' states
particles = initialStateRange(1) + (initialStateRange(2) - initialStateRange(1)) * rand(1, numParticles);
particles_history = zeros(T, numParticles);
particles_history(1, :) = particles;
particles_history_resamp = zeros(T, numParticles); % Storing the computed states for each particle accounting for their resamplings. The other history variable stores the histories that are compatible with the significant weights that got resampled for getting fed into the model.
particles_history_resamp(1, :) = particles;

particles_w = initialStateRange_w(1) + (initialStateRange_w(2) - initialStateRange_w(1)) * rand(1, numParticles);
particles_history_w = zeros(T, numParticles);
particles_history_w(1, :) = particles_w;
particles_history_resamp_w = zeros(T, numParticles); % Storing the computed states for each particle accounting for their resamplings. The other history variable stores the histories that are compatible with the significant weights that got resampled for getting fed into the model.
particles_history_resamp_w(1, :) = particles_w;

% Initialize all particles to start in mode one, reflecting 100% certainty
mode = ones(1, numParticles);  % Mode one for all particles
modeHistory = zeros(T, numParticles);
modeHistory(1, :) = mode; % Store the modes for the current timestep
modeHistory_resamp = zeros(T, numParticles); % Storing the modes for each particle accounting for their resamplings. The other history variable stores the histories that are compatible with the significant weights that got resampled for getting fed into the model.
modeHistory_resamp(1, :) = mode; % Store the modes for the current timestep

% Storing the final estimated states of the particle filter:
estimatedState_history = zeros(1,T);
estimatedState_history(1) = 0;

estimatedState_history_w = zeros(1,T);
estimatedState_history_w(1) = 0;

N_eff_history = zeros(1,T);
N_eff_history(1, 1) = numParticles;
N_eff_history_w = zeros(1,T);
N_eff_history_w(1, 1) = numParticles;

% Initialize weights uniformly
weights = ones(1, numParticles) / numParticles;
weights_history = zeros(T, numParticles);
weights_history(1, :) = weights;

weights_history_plt = zeros(T, numParticles);
weights_history_plt(1, :) = weights;

weights_w = ones(1, numParticles) / numParticles;
weights_history_w = zeros(T, numParticles);
weights_history_w(1, :) = weights_w;

weights_history_w_plt = zeros(T, numParticles);
weights_history_w_plt(1, :) = weights_w;


% Main filtering loop
for t = 2:1:T  % Assume T is the number of time steps
    measurement = getMeasurementAtTime(t, Z_real);  % Function to get measurement
    measurement_w = getMeasurementAtTime(t, Z_real_w);
    prior_output_all = zeros(2, numParticles);
    % Prediction step
    parfor p = 1:numParticles
        currentState = particles(p);
        currentMode = mode(p);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the importance density function for drawing the particles(p) from
        % or sampling them:

        % Parameters for the logistic distribution
        [prior1, prior2] = model_output(t, mode(p), particles_history, particles_history_w, modeHistory, p);
        prior_output = [prior1; prior2];


        particles(p) = (1-0.5)*(prior_output(1,:) + random(name, A_process_1, 1*B_process_1)) + 0.5*(measurement - random(name, A_measurement_1, 1*B_measurement_1));
        particles_w(p) = (1-0.5)*(prior_output(2,:) + random(name, A_process_1, 1*10*B_process_1)) + 0.5*(measurement_w - random(name, A_measurement_1, 1*10*B_measurement_1));

        prior_output_all(:, p) = prior_output;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Mode transition
        currentTransitions = transitionMatrix(currentMode, :);
        mode(p) = find(mnrnd(1, currentTransitions), 1);  % Sample new mode
    end


    particles_history(t, :) = particles;
    particles_history_resamp(t, :) = particles;
    particles_history_w(t, :) = particles_w;
    particles_history_resamp_w(t, :) = particles_w;
    modeHistory(t, :) = mode; % Store the modes for the current timestep
    modeHistory_resamp(t, :) = mode;



    % Measurement update (simplified for one state)
    
    weights = zeros(1, numParticles);
    parfor p = 1:numParticles

        % the importance density function for the weight update stage

        prior_output = prior_output_all(1, p);
        prior_output_w = prior_output_all(2, p);

        importance_density = kde_likelihood(particles(p), (1-0.5)*prior_output(1,:) + 0.5*measurement, kde_pdf_function_q);
        importance_density_w = kde_likelihood(particles_w(p), (1-0.5)*prior_output_w(1,:) + 0.5*measurement_w, kde_pdf_function_q_w);


        likelihood = pdf(name, measurement, particles(p) + A_measurement_1, B_measurement_1);
        likelihood_w = pdf(name, measurement_w, particles_w(p) + A_measurement_1, 10*B_measurement_1);

        % the process or prior distribution
        prior_x_i_minus_1 = prior_output;
        prior_x_i_minus_1_w = prior_output_w;

        prior = kde_likelihood(particles_history(t, p), prior_x_i_minus_1(1,:), kde_pdf_function);
        prior_w = kde_likelihood(particles_history_w(t, p), prior_x_i_minus_1_w(1,:), kde_pdf_function_w);

        % weight update
        weights(p) = likelihood*prior/importance_density*weights_history(t-1, p);

        weights_w(p) = likelihood_w*prior_w/importance_density_w*weights_history_w(t-1, p);

    end

    weights = weights / sum(weights);  % Normalize weights
    weights_w = weights_w / sum(weights_w);  % Normalize weights
    

    weights_history_plt(t, :) = weights;
    weights_history_w_plt(t, :) = weights_w;

    % Optional resampling step can be added here
    N_eff = 1/(sum(weights.^2, 'all'));
    N_eff_w = 1/(sum(weights_w.^2, 'all'));
    if N_eff > N_eff_thresh
        continue
    else
        [particles, mode, weights, particles_history, modeHistory] = resample(particles, mode, weights, numParticles, particles_history, modeHistory);
        [particles_w, weights_w, particles_history_w] = resample_w(particles_w, weights_w, numParticles, particles_history_w);
    end

    weights_history(t, :) = weights;
    weights_history_w(t, :) = weights_w;

    % Estimate state
    estimatedState = sum(particles .* weights);
    estimatedState_history(t) = estimatedState;

    estimatedState_w = sum(particles_w .* weights_w);
    estimatedState_history_w(t) = estimatedState_w;

    N_eff_history(1, t) = N_eff;
    N_eff_history_w(1, t) = N_eff_w;


    % Continue to next time step...
     %fprintf('Timestep: %d\n', t);
end



Z_filtered_markov = estimatedState_history;
end
%% Custom functions:

function measurement = getMeasurementAtTime(t, measurements)
% Check if t is within the range of measurements
if t >= 1 && t <= length(measurements)
    measurement = measurements(t);
else
    error('Time t is out of the range of available measurements');
end
end

function [x, m, w, p_h, m_h]=resample(x, m, w, N, p_h, m_h)
% Multinomial sampling with Ripley's method
u=cumprod(rand(1,N).^(1./[N:-1:1]));
u=fliplr(u);
wc=cumsum(w);
k=1;
for i=1:N
    while(wc(k)<u(i))
        k=k+1;
    end
    ind(i)=k;
end
x=x(:,ind);
m=m(:,ind);
p_h=p_h(:,ind);
%p_h_w = p_h_w(:,ind);
m_h=m_h(:,ind);
w=ones(1,N)./N;
%w_w=ones(1,N)./N;
end

function [x, w, p_h]=resample_w(x, w, N, p_h)
% Multinomial sampling with Ripley's method
u=cumprod(rand(1,N).^(1./[N:-1:1]));
u=fliplr(u);
wc=cumsum(w);
k=1;
for i=1:N
    while(wc(k)<u(i))
        k=k+1;
    end
    ind(i)=k;
end
x=x(:,ind);
p_h=p_h(:,ind);
%p_h_w = p_h_w(:,ind);
w=ones(1,N)./N;
%w_w=ones(1,N)./N;
end
