clear;
input_path = '../../inputs/inputs-1x-beta.xlsx';

inputs = readtable(input_path, 'Sheet', 'Pathogens');
inputs.gamma = 1 ./ inputs.infectious_period;
inputs.delta = 1 - ((1 - inputs.cfr) .^ inputs.gamma);

parameters = readtable(input_path, 'Sheet', 'Parameters');
population = parameters.population;
duration = parameters.duration;
replicates = 2;

initial_conditions = readtable(input_path, 'Sheet', 'Initial Conditions');

num_pathogens = height(inputs);

outputs = cell(replicates, 1);

% Run the simulation with a certain number of replicates
for replicate = 1:replicates
    % Start everyone out alive and allocate storage for each person for each
    % of the 4 compartments (S, E, I, R) per pathogen. Similar compartments
    % across pathogens are all grouped together (all S in one block of rows,
    % all E in the next block, etc.) Data types for each compartment are as
    % follows:
    %   - S: Bolean (0 or 1) indicated whether the individual is susceptible
    %   - E: Number of days the individual has been exposed
    %   - I: Number of days the individual has been infected
    %   - R: Number of days the individual has been recovered
    alive = true(1, population);
    compartments = zeros(num_pathogens * 4, population);
    
    % Declare variables for short-hand access to the various compartments
    S = (1 + 0 * (num_pathogens)):(1 * num_pathogens);
    E = (1 + 1 * (num_pathogens)):(2 * num_pathogens);
    I = (1 + 2 * (num_pathogens)):(3 * num_pathogens);
    R = (1 + 3 * (num_pathogens)):(4 * num_pathogens);
    
    % Specifify initial conditions for each compartment. Note that assignments
    % are made randomly and independently for each pathogen across the total
    % population.
    % compartments(I, 1:population) = rand(numel(I), population) < 0.05;
    % compartments(I, 1:population) = (compartments(I, 1:population) .* rand(numel(I), population)) .* inputs.infectious_period;
    % compartments(S, :) = ~logical(compartments(I, :));
    
    compartments = setInitialConditions(compartments, inputs, initial_conditions);
    
    % Allocate space for output data. The first group stores independent data
    % per compartment for each pathogen. The second set stores aggregate data
    % for total healthy, infected, and deceased individuals.
    outputs{replicate}.susceptible = zeros(num_pathogens, duration + 1);
    outputs{replicate}.exposures = zeros(num_pathogens, duration + 1);
    outputs{replicate}.infections = zeros(num_pathogens, duration + 1);
    outputs{replicate}.recoveries = zeros(num_pathogens, duration + 1);
    outputs{replicate}.fatalities = zeros(num_pathogens, duration + 1);
    
    outputs{replicate}.healthy = zeros(1, duration + 1);
    outputs{replicate}.infected = zeros(1, duration + 1);
    outputs{replicate}.deaths = zeros(1, duration + 1);
    
    % Add starting values to output (day = 0)
    susceptible = bitand(alive, compartments(S, :));
    exposed = compartments(E, :);
    infected = compartments(I, :);
    recovered = compartments(R, :);
    
    outputs{replicate}.susceptible(:, 1) = sum(susceptible > 0, 2);
    outputs{replicate}.exposures(:, 1) = sum(exposed > 0, 2);
    outputs{replicate}.infections(:, 1) = sum(infected > 0, 2);
    outputs{replicate}.recoveries(:, 1) = sum(recovered > 0, 2);
    outputs{replicate}.fatalities(:, 1) = 0;
    
    outputs{replicate}.healthy(1) = sum(sum(susceptible, 1) == num_pathogens);
    outputs{replicate}.infected(1) = sum(logical(sum(infected > 0, 1)));
    outputs{replicate}.deaths(1) = sum(~alive);
    
    for day=2:duration + 1
        N = sum(alive);
        
        % Deconstruct intial values for each compartment. The susceptible
        % population is masked by the alive flag to ignore any susceptibility
        % states for deceased persons.
        susceptible = bitand(alive, compartments(S, :));
        exposed = compartments(E, :);
        infected = compartments(I, :);
        recovered = compartments(R, :);
        
        % Determine the instantaneous lambda value based on the current number
        % of infections for each pathogen.
        infections = bitand(alive, compartments(I, :) > 0);
        total_infections = sum(infections, 2);
        lambda = (inputs.beta .* total_infections) / N;
        
        % Draw random variates to determine which individuals will become
        % exposed per pathogen. Exposed indviduals will also have their
        % susceptible flag set to false.
        variates = rand(num_pathogens, population);
        newly_exposed = susceptible & (variates < lambda);
        susceptible(newly_exposed) = false;
        exposed = logical(exposed) .* (exposed + 1);
        exposed(newly_exposed) = 1;
        
        % Identify individuals that should be moved from E to I based on the
        % passage of the appropriate latency for each pathogen (in days).
        % Persons that move to infected, will have their E reset to 0 days.
        newly_infected = exposed > inputs.latency;
        exposed(newly_infected) = 0;
        infected = bitand(alive, logical(infected)) .* (infected + 1);
        infected(newly_infected) = 1;
        
        % Draw random variates to determine identify any disease-related
        % fatalities per pathogen. It's possible for a person to die
        % simultaneously from a co-infection. This is considereld a single
        % death and then used to update the alive flag for those individuals.
        % The susceptible flag associated with the pathogen responsible for
        % death is also set to true to track specific pathogen-related deaths.
        variates = rand(num_pathogens, population);
        fatalities = logical(infected) & (variates < inputs.delta);
        deaths = sum(fatalities, 1) > 0;
        alive(deaths) = false;
        susceptible(fatalities) = true;
        
        % Identify individuals that should be moved from I to R based on the
        % passage of the appropriate infectious period for each pathogen (in
        % days). Persons that move to recovered, will have their I reset to
        % 0 days.
        newly_recovered = infected > inputs.infectious_period;
        infected(newly_recovered) = 0;
        recovered = logical(recovered) .* (recovered + 1);
        recovered(newly_recovered) = 1;
        
        % Identify individuals that should be moved from R to S based on the
        % passage of the appropriate immunity period for each pathogen (in
        % days). Persons that move to susceptible, will have their R reset to
        % 0 days.
        newly_susceptible = recovered > inputs.immunity_period;
        susceptible(newly_susceptible) = true;
        recovered(newly_susceptible) = 0;
        
        % Update the compartments variable with the changes from this
        % iteration, as well as the output data.
        compartments(S, :) = susceptible;
        compartments(E, :) = exposed;
        compartments(I, :) = infected;
        compartments(R, :) = recovered;
        
        outputs{replicate}.susceptible(:, day) = sum(susceptible > 0, 2);
        outputs{replicate}.exposures(:, day) = sum(exposed > 0, 2);
        outputs{replicate}.infections(:, day) = sum(infected > 0, 2);
        outputs{replicate}.recoveries(:, day) = sum(recovered > 0, 2);
        outputs{replicate}.fatalities(:, day) = outputs{replicate}.fatalities(:, max(1, day - 1)) + sum(fatalities, 2);
        
        outputs{replicate}.healthy(day) = sum(sum(susceptible, 1) == num_pathogens);
        outputs{replicate}.infected(day) = sum(logical(sum(infected > 0, 1)));
        outputs{replicate}.deaths(day) = sum(~alive);
    end
    
    disp(['Complete run ', num2str(replicate), '...']);
end

days = 0:duration;

outputs_mat = cell2mat(outputs);


figure('Name', 'Agent SEIRDS - All');
for pathogen = 1:num_pathogens
    subplot(ceil(num_pathogens / 3), 3, pathogen);
    plot(days, [ ...
        outputs{1}.susceptible(pathogen, :); ...
        outputs{1}.exposures(pathogen, :); ...
        outputs{1}.infections(pathogen, :); ...
        outputs{1}.recoveries(pathogen, :); ...
        outputs{1}.fatalities(pathogen, :)], ...
        'LineWidth', 2);
    title(inputs.pathogen(pathogen));
    legend('S', 'E', 'I', 'R', 'D');
end

figure('Name', 'Agent SEIRDS - Totals');
plot(days, [ ...
    outputs{1}.healthy; ...
    outputs{1}.infected; ...
    outputs{1}.deaths], ...
    'LineWidth', 2);
legend('Healthy', 'Infected', 'Deaths');
title("Total diarrhea");

function initials = setInitialConditions(compartments, inputs, conditions)
%% Sets initial conditions for the compartment model
%
% Inputs:
%   - compartments: The initial matrix of compartment data or all
%       pathogens.

num_pathogens = height(conditions);
population = size(compartments, 2);

total = sum([conditions.e, conditions.i, conditions.r], 'all');
if (total > 1)
    ME = MException('setInitialConditions:inputError','Initial conditions must sum to less than unit.');
    throw(ME);
end

S = (1 + 0 * (num_pathogens)):(1 * num_pathogens);
E = (1 + 1 * (num_pathogens)):(2 * num_pathogens);
I = (1 + 2 * (num_pathogens)):(3 * num_pathogens);
R = (1 + 3 * (num_pathogens)):(4 * num_pathogens);

person = 1;
for pathogen = 1:num_pathogens
    condition = conditions(pathogen, :);
    input = inputs(pathogen, :);
    initial_person = person;
    
    exposedPersons = round(condition.e * population);
    infectedPersons = round(condition.i * population);
    recoveredPersons = round(condition.r * population);
    
    if (exposedPersons > 0)
        range = person:(person + exposedPersons - 1);
        compartments(E(1) + pathogen - 1, range) = ...
            round(rand(1, numel(range)) .* input.latency);
        person = range(end) + 1;
    end
    
    if (infectedPersons > 0)
        range = person:(person + infectedPersons - 1);
        compartments(I(1) + pathogen - 1, range) = ...
            round(rand(1, numel(range)) .* input.infectious_period);
        person = range(end) + 1;
    end
    
    if (recoveredPersons > 0)
        range = person:(person + recoveredPersons - 1);
        compartments(R(1) + pathogen - 1, range) = ...
            round(rand(1, numel(range)) .* input.immunity_period);
        person = range(end) + 1;
    end
    
    compartments(S, initial_person:person) = 1;
    compartments(S(1) + pathogen - 1, initial_person:person) = 0;
end

if (person <= population)
    compartments(S, person:population) = 1;
end

initials = compartments;
end