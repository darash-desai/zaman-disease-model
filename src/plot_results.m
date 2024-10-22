%% Plots results from a recently run model.

% The path to the input file used to run the model
input_path = 'inputs-baseline.xlsx';

inputs = readtable(input_path, 'Sheet', 'Pathogens');
parameters = readtable(input_path, 'Sheet', 'Parameters');

days = 0:parameters.duration;
num_pathogens = height(inputs);

population = parameters.population;

figure(1);
for pathogen = 1:num_pathogens
    data = pathogen_data(pathogen);
    susceptible = data.susceptible;
    exposures = data.exposures;
    infections = data.infections;
    recoveries = data.recoveries;
    % new_infections = data.new_infections;
    fatalities = data.fatalities;
    
    subplot(ceil(num_pathogens / 3), 3, pathogen);
    hold on;
    ploterrors(days, 100 * susceptible / population, '#77AC30');
    ploterrors(days, 100 * exposures / population, '#EDB120');
    ploterrors(days, 100 * infections / population, '#D95319');
    ploterrors(days, 100 * recoveries / population, '#0072BD');
    % ploterrors(days, 100 * cumsum(data.new_infections, 2) / population, '#7E2F8E');
    ploterrors(days, 100 * cumsum(fatalities, 2) / population, '#A2142F');
    
    title(inputs.pathogen(pathogen));
    xlabel('Days');
    ylabel('Population (%)');
    % legend('S', 'E', 'I', 'R', 'Cases', 'D');
    legend('S', 'E', 'I', 'R', 'D');
end

figure(2);
hold on;

healthy = summary_data.healthy;
infected = summary_data.infected;
new_cases = summary_data.new_cases;
deaths = summary_data.deaths;

ploterrors(days, 100 * healthy / population, '#77AC30');
ploterrors(days, 100 * infected / population, '#D95319');
ploterrors(days, 100 * cumsum(deaths / population, 2), '#A2142F');

legend('Healthy', 'Infected', 'Deaths (cumulative)');
title("Total diarrhea");
xlabel('Days');
ylabel('Population (%)');

figure(3);
hold on;
ploterrors(days, 100 * cumsum(new_cases, 2) / population, '#7E2F8E');
title("Total diarrhea cases");
xlabel('Days');
ylabel('Population (%)');