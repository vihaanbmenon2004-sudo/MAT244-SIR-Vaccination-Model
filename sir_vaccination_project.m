%% MAT244 Final Project
% Extended SIR Model with Vaccination
% Vihaan Menon, Crishna Divekar, Farnam Hesabi

clear;
clc;
close all;

%% ============================================================
% 1. BASELINE PARAMETERS
% =============================================================

gamma = 1/7;             % recovery rate (1/day)
R0 = 2.5;                % basic reproduction number
beta = R0 * gamma;       % transmission rate

S0 = 0.999;              % initially susceptible
I0 = 0.001;              % initially infected
Rinit = 0;               % initially recovered
C0 = I0;                 % cumulative infections

y0 = [S0; I0; Rinit; C0];

tspan = [0 180];

options = odeset('RelTol',1e-9,'AbsTol',1e-11);

fprintf('Baseline parameters:\n');
fprintf('R0 = %.2f\n',R0);
fprintf('beta = %.6f per day\n',beta);
fprintf('gamma = %.6f per day\n\n',gamma);


%% ============================================================
% 2. BASELINE SIR MODEL: NO VACCINATION
% =============================================================

nu0 = 0;
delay = 0;

[t,y] = ode45(@(t,y) sirVaccinationODE(t,y,beta,gamma,nu0,delay), ...
              tspan,y0,options);

S = y(:,1);
I = y(:,2);
R = y(:,3);
C = y(:,4);

[Imax,indexPeak] = max(I);
tPeak = t(indexPeak);
Cfinal = C(end);

epsilon = 1e-5;

afterPeak = find(t > tPeak & I < epsilon,1);

if isempty(afterPeak)
    epidemicDuration = NaN;
else
    epidemicDuration = t(afterPeak);
end

conservationError = max(abs(S + I + R - 1));

fprintf('BASELINE RESULTS\n');
fprintf('Peak infected = %.4f (%.2f%%)\n',Imax,100*Imax);
fprintf('Peak day = %.2f\n',tPeak);
fprintf('Cumulative infected = %.4f (%.2f%%)\n',Cfinal,100*Cfinal);
fprintf('Epidemic duration = %.2f days\n',epidemicDuration);
fprintf('Maximum conservation error = %.3e\n\n',conservationError);


%% ============================================================
% 3. IMMEDIATE VACCINATION RATE COMPARISON
% =============================================================

vaccinationRates = [0 0.005 0.010 0.020 0.030 0.050];

peakInfected = zeros(size(vaccinationRates));
peakDay = zeros(size(vaccinationRates));
cumulativeInfected = zeros(size(vaccinationRates));

figure;
hold on;

for k = 1:length(vaccinationRates)

    nu0 = vaccinationRates(k);
    delay = 0;

    [t,y] = ode45(@(t,y) sirVaccinationODE( ...
        t,y,beta,gamma,nu0,delay),tspan,y0,options);

    I = y(:,2);
    C = y(:,4);

    [peakInfected(k),idx] = max(I);
    peakDay(k) = t(idx);
    cumulativeInfected(k) = C(end);

    % Plot only rates used in main figure
    if nu0 <= 0.03
        plot(t,I,'LineWidth',1.5, ...
            'DisplayName',sprintf('\\nu = %.3f',nu0));
    end
end

xlabel('Time (days)');
ylabel('Infected proportion');
title('Effect of Immediate Vaccination Rate on Infection Prevalence');
legend('Location','best');
grid on;
hold off;

saveas(gcf,'figure1_vaccination_rates.png');


fprintf('IMMEDIATE VACCINATION RESULTS\n');

fprintf('%8s %15s %12s %20s\n', ...
    'nu','Peak infected','Peak day','Cumulative infected');

for k = 1:length(vaccinationRates)

    fprintf('%8.3f %14.2f%% %12.2f %19.2f%%\n', ...
        vaccinationRates(k), ...
        100*peakInfected(k), ...
        peakDay(k), ...
        100*cumulativeInfected(k));

end

fprintf('\n');


%% ============================================================
% 4. DELAYED VACCINATION
% =============================================================

nu0 = 0.020;

delays = [0 10 20 30 40];

delayPeak = zeros(size(delays));
delayPeakDay = zeros(size(delays));
delayCumulative = zeros(size(delays));

figure;
hold on;

for k = 1:length(delays)

    delay = delays(k);

    [t,y] = ode45(@(t,y) sirVaccinationODE( ...
        t,y,beta,gamma,nu0,delay),tspan,y0,options);

    I = y(:,2);
    C = y(:,4);

    [delayPeak(k),idx] = max(I);
    delayPeakDay(k) = t(idx);
    delayCumulative(k) = C(end);

    plot(t,I,'LineWidth',1.5, ...
        'DisplayName',sprintf('Delay = %d days',delay));

end

xlabel('Time (days)');
ylabel('Infected proportion');
title('Effect of Vaccination Delay on Infection Prevalence');
legend('Location','best');
grid on;
hold off;

saveas(gcf,'figure2_vaccination_delays.png');


fprintf('VACCINATION DELAY RESULTS\n');

fprintf('%10s %15s %12s %20s\n', ...
    'Delay','Peak infected','Peak day','Cumulative infected');

for k = 1:length(delays)

    fprintf('%8d %14.2f%% %12.2f %19.2f%%\n', ...
        delays(k), ...
        100*delayPeak(k), ...
        delayPeakDay(k), ...
        100*delayCumulative(k));

end

fprintf('\n');


%% ============================================================
% 5. TWO-PARAMETER SWEEP: beta vs vaccination rate
% =============================================================

betaValues = linspace(0.20,0.50,31);
nuValues = linspace(0,0.04,31);

peakMatrix = zeros(length(nuValues),length(betaValues));

for i = 1:length(nuValues)

    for j = 1:length(betaValues)

        betaTest = betaValues(j);
        nuTest = nuValues(i);

        [~,y] = ode45(@(t,y) sirVaccinationODE( ...
            t,y,betaTest,gamma,nuTest,0), ...
            tspan,y0,options);

        I = y(:,2);

        peakMatrix(i,j) = max(I);

    end

end

figure;

imagesc(betaValues,nuValues,peakMatrix);

set(gca,'YDir','normal');

colorbar;

xlabel('Transmission rate \beta');
ylabel('Vaccination rate \nu');
title('Peak Infection Proportion Across Transmission and Vaccination Rates');

saveas(gcf,'figure3_heatmap.png');


%% ============================================================
% 6. ±20% PARAMETER SENSITIVITY ANALYSIS
% =============================================================

baselineNu = 0.020;

[tBase,yBase] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,beta,gamma,baselineNu,0),tspan,y0,options);

baselinePeak = max(yBase(:,2));
baselineCumulative = yBase(end,4);

labels = { ...
    '\beta -20%', ...
    '\beta +20%', ...
    '\gamma -20%', ...
    '\gamma +20%', ...
    '\nu -20%', ...
    '\nu +20%'};

peakSensitivity = zeros(6,1);
cumSensitivity = zeros(6,1);

%% beta -20%
betaTest = 0.8*beta;

[~,y] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,betaTest,gamma,baselineNu,0),tspan,y0,options);

peakSensitivity(1) = max(y(:,2));
cumSensitivity(1) = y(end,4);


%% beta +20%
betaTest = 1.2*beta;

[~,y] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,betaTest,gamma,baselineNu,0),tspan,y0,options);

peakSensitivity(2) = max(y(:,2));
cumSensitivity(2) = y(end,4);


%% gamma -20%
gammaTest = 0.8*gamma;

[~,y] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,beta,gammaTest,baselineNu,0),tspan,y0,options);

peakSensitivity(3) = max(y(:,2));
cumSensitivity(3) = y(end,4);


%% gamma +20%
gammaTest = 1.2*gamma;

[~,y] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,beta,gammaTest,baselineNu,0),tspan,y0,options);

peakSensitivity(4) = max(y(:,2));
cumSensitivity(4) = y(end,4);


%% vaccination -20%
nuTest = 0.8*baselineNu;

[~,y] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,beta,gamma,nuTest,0),tspan,y0,options);

peakSensitivity(5) = max(y(:,2));
cumSensitivity(5) = y(end,4);


%% vaccination +20%
nuTest = 1.2*baselineNu;

[~,y] = ode45(@(t,y) sirVaccinationODE( ...
    t,y,beta,gamma,nuTest,0),tspan,y0,options);

peakSensitivity(6) = max(y(:,2));
cumSensitivity(6) = y(end,4);


%% Calculate percentage changes

peakPercentChange = ...
    100*(peakSensitivity-baselinePeak)/baselinePeak;

cumPercentChange = ...
    100*(cumSensitivity-baselineCumulative)/baselineCumulative;


fprintf('SENSITIVITY ANALYSIS\n');

fprintf('%15s %15s %15s %18s %18s\n', ...
    'Change','Peak','Peak change','Cumulative','Cum. change');

for k = 1:6

    fprintf('%15s %14.2f%% %14.2f%% %17.2f%% %17.2f%%\n', ...
        labels{k}, ...
        100*peakSensitivity(k), ...
        peakPercentChange(k), ...
        100*cumSensitivity(k), ...
        cumPercentChange(k));

end


%% Sensitivity plot

figure;

bar(peakPercentChange);

yline(0,'k-');

set(gca,'XTick',1:6);
set(gca,'XTickLabel',labels);
xtickangle(25);

ylabel('Change in peak infection (%)');
title('One-at-a-Time Sensitivity of Peak Infection');
grid on;

saveas(gcf,'figure4_sensitivity.png');


%% ============================================================
% LOCAL FUNCTION: EXTENDED SIR MODEL
% =============================================================

function dydt = sirVaccinationODE(t,y,beta,gamma,nu0,delay)

    S = y(1);
    I = y(2);

    % Vaccination begins after specified delay
    if t >= delay
        nu = nu0;
    else
        nu = 0;
    end

    infection = beta*S*I;

    dSdt = -infection - nu*S;
    dIdt = infection - gamma*I;
    dRdt = gamma*I + nu*S;

    % cumulative number ever infected
    dCdt = infection;

    dydt = [dSdt; dIdt; dRdt; dCdt];

end