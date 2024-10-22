%% Traditional implementation of an SEIRDS model

function outputs = seirds( ...
    duration, ...
    beta, ...
    latency, ...
    infectious_period, ...
    CFR, ...
    immunity_period, ...
    population, ...
    initial_proportions, ...
    birth_rate, ...
    death_rate ...
    )
%% SEIRDS model
% Calculates the number of individuals that mirate through the various
% compartments of the SEIRDS model over the duration of days specified.
%
% Parameters:
%   - duration: The duration of the model, in days.
%   - beta: The transmission probability when exposed to an infectious
%       person.
%   - latency: The average latency period, in days.
%   - infectious_period: The average infectious period, in days.
%   - CFR: The case fatality rate.
%   - immunity_period: The average duration of immunity, in days.
%   - population: The total population to consider.
%   - initial_proportions: The initial proportion of the total
%       population in each compartment (s, e, i, r, d).

% The infectious rate, describing the rate at which exposed individuals
% become infectous. The reciprocal of this rate describes the average
% latency period in days.
sigma = 1 / latency;

% The recovery rate, describing the rate at which infectious individuals
% recover. The reciprocal of this rate describes the average duration of
% infectiousness in days.
gamma = 1 / infectious_period;

% The disease-related death rate, describing probability of an individual
% dying on a given day based on the case fatality rate and average duration
% of infectiousness.
delta = 1 - ((1 - CFR) ^ gamma);

% The immunity wane rate, describing the rate at which recovered
% individuals lose their immunity and become susceptible again. The
% reciprocal of this is the average duration it takes to lose immunity in
% days.
eta = 1 / immunity_period;

% The regional natural birth rate in persons per day.
mu = birth_rate;

% The regional natural death rate in persons per day.
nu = death_rate;

% Check to make sure the initial proportions sum to unity.
if sum(initial_proportions) ~= 1
    disp("Error: Invalid initial conditions.");
    return;
end

compartments = num2cell(population * initial_proportions);
[S, E, I, R, D] = compartments{:};

outputs.s = zeros([duration + 1, 1]);
outputs.e = zeros([duration + 1, 1]);
outputs.i = zeros([duration + 1, 1]);
outputs.r = zeros([duration + 1, 1]);
outputs.d = zeros([duration + 1, 1]);
outputs.ni = zeros([duration + 1, 1]);

outputs.s(1) = S;
outputs.e(1) = E;
outputs.i(1) = I;
outputs.r(1) = R;
outputs.d(1) = D;
outputs.ni(1) = 0;

% Function defining the probability of being exposed to an infectious
% person in the total population.
exposure_prob = @(I, N) I/N;

for day = 2:duration + 1
    N = S + E + I + R;
    lambda = beta * exposure_prob(I, N);
    new_infections = sigma * E;
    
    dD = delta * I;
    dS = mu * N + eta * R - (lambda  + nu) * S;
    dE = lambda * S - sigma * E - nu * E;
    dI = new_infections - (gamma + nu) * I - dD;
    dR = gamma * I - eta * R - nu * R;
    
    S = max(0, S + dS);
    E = max(0, E + dE);
    I = max(0, I + dI);
    R = max(0, R + dR);
    D = D + dD;
    
    outputs.s(day) = S;
    outputs.e(day) = E;
    outputs.i(day) = I;
    outputs.r(day) = R;
    outputs.d(day) = D;
    outputs.ni(day) = new_infections;
end
