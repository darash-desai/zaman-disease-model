% Resource: https://goodcalculators.com/seir-model-calculator/

clear;
input_path = 'inputs-baseline-density.xlsx';

inputs = readtable(input_path, 'Sheet', 'Pathogens');
parameters = readtable(input_path, 'Sheet', 'Parameters');
initial_conditions = readtable(input_path, 'Sheet', 'Initial Conditions');

inputs.gamma = 1 ./ inputs.infectious_period;
inputs.delta = 1 - ((1 - inputs.cfr) .^ inputs.gamma);

% Adjust infectious period/gamma
inputs.infectious_period = round(inputs.infectious_period * 1.5);

population = parameters.population;
duration = parameters.duration;
replicates = parameters.replicates;
low_density = parameters.low_density;
med_density = parameters.med_density;
high_density = parameters.high_density;
low_mixing = parameters.low_mixing;
med_mixing = parameters.med_mixing;
high_mixing = parameters.high_mixing;
i_stdev = parameters.i_stdev;
r_stdev = parameters.r_stdev;
lambda_exp = parameters.lambda_exp;

num_pathogens = height(inputs);

blank_struct = @(A) struct('susceptible', zeros(replicates, duration + 1), ...
                           'exposures', zeros(replicates, duration + 1), ...
                           'infections', zeros(replicates, duration + 1), ...
                           'new_infections', zeros(replicates, duration + 1), ...
                           'recoveries', zeros(replicates, duration + 1), ...
                           'fatalities', zeros(replicates, duration + 1));
pathogen_data = cellfun(blank_struct, cell(1,num_pathogens));

summary_data = struct;
summary_data.healthy = zeros(replicates, duration + 1);
summary_data.infected = zeros(replicates, duration + 1);
summary_data.new_cases = zeros(replicates, duration + 1);
summary_data.deaths = zeros(replicates, duration + 1);

% Run the simulation with a certain number of replicates
for replicate = 1:replicates    
    % Start everyone out alive and allocate storage for each person for each 
    % of the 4 compartments (S, E, I, R), per pathogen. Similar compartments 
    % across pathogens are all grouped together (all S in one block of rows, 
    % all E in the next block, etc.) Data types for each compartment are as 
    % follows:
    %   - S: Bolean (0 or 1) indicated whether the individual is susceptible
    %   - E: Number of days the individual has been exposed
    %   - I: Number of days the individual has been infected
    %   - R: Number of days the individual has been recovered
    alive = true(1, population);
    mixing_param = ones(1, population);
    compartments = zeros(num_pathogens * 4, population);

    % Allocate space for personalized infectious periods
    infectious_periods = zeros(num_pathogens, population);
    for pathogen = 1:num_pathogens
        mu = inputs.infectious_period(pathogen);
        sigma = mu * i_stdev;        
        infectious_periods(pathogen, :) = ...
            max(1, round(sigma .* randn(1, population) + mu));
    end

    % Allocate space for personalized recovery periods
    recovery_periods = zeros(num_pathogens, population);
    for pathogen = 1:num_pathogens
        mu = inputs.immunity_period(pathogen);
        sigma = mu * r_stdev;        
        recovery_periods(pathogen, :) = ...
            max(1, round(sigma .* randn(1, population) + mu));
    end
    
    % Declare variables for short-hand access to the various compartments
    S = (1 + 0 * (num_pathogens)):(1 * num_pathogens);
    E = (1 + 1 * (num_pathogens)):(2 * num_pathogens);
    I = (1 + 2 * (num_pathogens)):(3 * num_pathogens);
    R = (1 + 3 * (num_pathogens)):(4 * num_pathogens);

    % Declare variables for short-hand access to density groups
    LD = 1:round(population * low_density);
    MD = LD(end):LD(end) + round(population * med_density);
    HD = MD(end):population;
    
    mixing_param(LD) = low_mixing;
    mixing_param(MD) = med_mixing;
    mixing_param(HD) = high_mixing;
    mixing_param = mixing_param(randperm(length(mixing_param)));

    % Specifify initial conditions for each compartment. Note that assignments
    % are made randomly and independently for each pathogen across the total
    % population.
    initials = setInitialConditions(compartments, inputs, initial_conditions, infectious_periods, recovery_periods);
    compartments = initials;
    
    % Allocate space for output data. The first group stores independent data
    % per compartment for each pathogen. The second set stores aggregate data
    % for total healthy, infected, and deceased individuals.
    outputs_susceptible = zeros(num_pathogens, duration + 1);
    outputs_exposures = zeros(num_pathogens, duration + 1);
    outputs_new_infections = zeros(num_pathogens, duration + 1);
    outputs_infections = zeros(num_pathogens, duration + 1);
    outputs_recoveries = zeros(num_pathogens, duration + 1);
    outputs_fatalities = zeros(num_pathogens, duration + 1);
    
    outputs_healthy = zeros(1, duration + 1);
    outputs_infected = zeros(1, duration + 1);
    outputs_new_cases = zeros(1, duration + 1);
    outputs_deaths = zeros(1, duration + 1);
    
    % Add starting values to output (day = 0)
    susceptible = bitand(alive, compartments(S, :));
    exposed = compartments(E, :);
    infected = compartments(I, :);
    recovered = compartments(R, :);
    
    outputs_susceptible(:, 1) = sum(susceptible > 0, 2);
    outputs_exposures(:, 1) = sum(exposed > 0, 2);
    outputs_infections(:, 1) = sum(infected > 0, 2);
    outputs_recoveries(:, 1) = sum(recovered > 0, 2);
    outputs_fatalities(:, 1) = 0;
    
    outputs_healthy(1) = sum(sum(susceptible, 1) == num_pathogens);
    outputs_infected(1) = sum(logical(sum(infected > 0, 1)));
    outputs_new_cases(1) = 0;
    outputs_deaths(1) = sum(~alive);
    
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
        lambda = (inputs.beta .* ((total_infections / N) .^ mixing_param)) .^ lambda_exp;
        
        % Draw random variates to determine which individuals will become
        % exposed per pathogen. Exposed indviduals will also have their
        % susceptible flag set to false.
        variates = rand(num_pathogens, population);
        infectable = ones(num_pathogens, population);
        infected_persons = logical(sum(infections, 1));
        
        % Adjust probability of infection if you're already infected with
        % a pathogen. This reduces lambda by this factor.
        infectable(:, infected_persons) = 1;

        newly_exposed = susceptible & (variates < infectable .* lambda);
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
    
        % Draw random variates to identify any disease-related
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
        newly_recovered = infected > infectious_periods;
        infected(newly_recovered) = 0;
        recovered = logical(recovered) .* (recovered + 1);
        recovered(newly_recovered) = 1;
    
        % Identify individuals that should be moved from R to S based on the
        % passage of the appropriate immunity period for each pathogen (in 
        % days). Persons that move to susceptible, will have their R reset to 
        % 0 days.
        newly_susceptible = recovered > recovery_periods;
        susceptible(newly_susceptible) = true;
        recovered(newly_susceptible) = 0;
    
        % Update the compartments variable with the changes from this 
        % iteration, as well as the output data.
        compartments(S, :) = susceptible;
        compartments(E, :) = exposed;
        compartments(I, :) = infected;
        compartments(R, :) = recovered;
        
        outputs_susceptible(:, day) = sum(susceptible > 0, 2);
        outputs_exposures(:, day) = sum(exposed > 0, 2);
        outputs_infections(:, day) = sum(infected > 0, 2);
        outputs_new_infections(:, day) = sum(newly_infected, 2);
        outputs_recoveries(:, day) = sum(recovered > 0, 2);
        outputs_fatalities(:, day) = sum(fatalities, 2);
    
        outputs_healthy(day) = sum(sum(susceptible, 1) == num_pathogens);
        outputs_infected(day) = sum(logical(sum(infected > 0, 1)));
        outputs_new_cases(day) = sum(sum(newly_infected, 1) > 0);
        outputs_deaths(day) = sum(deaths, 2);
    end

    % Reorganize output data by pathogen
    for pathogen = 1:num_pathogens
        pathogen_data(pathogen).susceptible(replicate, :) = outputs_susceptible(pathogen, :);
        pathogen_data(pathogen).exposures(replicate, :) = outputs_exposures(pathogen, :);
        pathogen_data(pathogen).infections(replicate, :) = outputs_infections(pathogen, :);
        pathogen_data(pathogen).new_infections(replicate, :) = outputs_new_infections(pathogen, :);
        pathogen_data(pathogen).recoveries(replicate, :) = outputs_recoveries(pathogen, :);
        pathogen_data(pathogen).fatalities(replicate, :) = outputs_fatalities(pathogen, :);
    end

    summary_data.healthy(replicate, :) = outputs_healthy;
    summary_data.infected(replicate, :) = outputs_infected;
    summary_data.new_cases(replicate, :) = outputs_new_cases;
    summary_data.deaths(replicate, :) = outputs_deaths;

    disp(['Completed run ', num2str(replicate), '...']);
end

days = 0:duration;

figure('Name', 'Agent SEIRDS - All');
for pathogen = 1:num_pathogens
    data = pathogen_data(pathogen);

    subplot(ceil(num_pathogens / 3), 3, pathogen);
    hold on;
    ploterrors(days, data.susceptible, '#77AC30');
    ploterrors(days, data.exposures, '#EDB120');
    ploterrors(days, data.infections, '#D95319');
    ploterrors(days, data.recoveries, '#0072BD');
    ploterrors(days, cumsum(data.new_infections, 2), '#7E2F8E');
    ploterrors(days, cumsum(data.fatalities, 2), '#A2142F');

    title(inputs.pathogen(pathogen));
    xlabel('Days');
    ylabel('Population');
    legend('S', 'E', 'I', 'R', 'Cases', 'D');
end

figure('Name', 'Agent SEIRDS - Totals');
hold on;

ploterrors(days, summary_data.healthy, '#77AC30');
ploterrors(days, summary_data.infected, '#D95319');
ploterrors(days, cumsum(summary_data.new_cases, 2), '#7E2F8E');
ploterrors(days, cumsum(summary_data.deaths, 2), '#A2142F');

legend('Healthy', 'Infected', 'Cases (cumlative)', 'Deaths (cumulative)');
title("Total diarrhea");
xlabel('Days');
ylabel('Population');

function initials = setInitialConditions(compartments, inputs, conditions, infectious_periods, recovery_periods)
    %% Sets initial conditions for the compartment model
    %
    % This function sets up the initial configuration for the E, I, and R
    % compartments of the SEIRS model given the proportion of agents that  
    % should be attributed to each compartment for each pathogen. It treats
    % each pathogen independently when assigning individuals to the E, I,
    % and R compartments, so co-infections are possible as are agents that
    % are in the immunity period for several pathogens.
    %
    % Inputs:
    %   - compartments: The initial matrix of compartment data or all
    %       pathogens.    
    %   - inputs: The pathogen-related inputs.
    %   - conditions: The initial conditions for the proportion of the
    %       population in each compartment broken out by pathogen.

    num_pathogens = height(conditions);
    population = size(compartments, 2);    

    S = (1 + 0 * (num_pathogens)):(1 * num_pathogens);
    E = (1 + 1 * (num_pathogens)):(2 * num_pathogens);
    I = (1 + 2 * (num_pathogens)):(3 * num_pathogens);
    R = (1 + 3 * (num_pathogens)):(4 * num_pathogens);

    % Start all persons as susceptible to all pathogens
    compartments(S, :) = 1;

    for pathogen = 1:num_pathogens
        condition = conditions(pathogen, :);
        input = inputs(pathogen, :);

        total = sum([condition.e, condition.i, condition.r]);
        if (total > 1)
            ME = MException('setInitialConditions:inputError','Initial conditions must sum to less than unit.');
            throw(ME);
        end

        exposedPersons = round(condition.e * population);
        infectedPersons = round(condition.i * population);
        recoveredPersons = round(condition.r * population);

        range = randperm(population, exposedPersons + infectedPersons + recoveredPersons);

        if (exposedPersons > 0)
            exposedRange = range(1:exposedPersons);
            compartments(E(1) + pathogen - 1, exposedRange) = ...
                randi(input.latency, 1, numel(exposedRange));

            % Exposed persons should no longer be susceptible
            compartments(S(1) + pathogen - 1, exposedRange) = 0;
        end

        if (infectedPersons > 0)
            infectedRange = range((exposedPersons + 1):(exposedPersons + infectedPersons));
            compartments(I(1) + pathogen - 1, infectedRange) = ...
                1 + round(rand(1, numel(infectedRange)) .* (infectious_periods(pathogen, infectedRange) - 1));

            % Infected persons should no longer be susceptible
            compartments(S(1) + pathogen - 1, infectedRange) = 0;
        end

        if (recoveredPersons > 0)
            recoveredRange = range((exposedPersons + infectedPersons + 1):numel(range));
            compartments(R(1) + pathogen - 1, recoveredRange) = ...
                1 + round(rand(1, numel(recoveredRange)) .* (recovery_periods(pathogen, recoveredRange) - 1));

            % Recovered persons should no longer be susceptible
            compartments(S(1) + pathogen - 1, recoveredRange) = 0;
        end

        % Check
        exposed = 100 * sum(logical(compartments(E(1) + pathogen - 1, :))) / population;
        infected = 100 * sum(logical(compartments(I(1) + pathogen - 1, :))) / population;
        recovered = 100 * sum(logical(compartments(R(1) + pathogen - 1, :))) / population;
        disp(['Stats for ' input.pathogen{1} ' [e: ' num2str(exposed) '%, i: ' num2str(infected) '%, r:' num2str(recovered) '%]']);
    end

    initials = compartments;
end



% TODO
%
% 1. 