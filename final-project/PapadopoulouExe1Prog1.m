% Analush kai Epeksergasia Dedomenwn
% Zhthma 1 
% Papadopoulou Martha
% AEM: 4438

clear; clc; close all;

% import data
data = importdata('SeoulBike.xlsx');
bikes = data.data(:,1); %rented bike count
season = data.data(:,end-1);

% plot histograms of number of bikes rented for each season
figure(1)
season_names = {'Winter', 'Spring', 'Summer', 'Autumn'};

for s = 1:4
    % access the right season's bike data
    season_bikes = bikes(season == s);
    
    % create subplots
    subplot(2, 2, s);
    hold on;
    
    % plot histogram
    histogram(season_bikes);
    xlabel('Number of Bikes Rented');
    ylabel('Frequency');
    title(['Season ', season_names{s}]);
end

%%
%try different distributions
figure(2);

for s = 1:4
    % access the right season's bike data
    season_bikes = bikes(season == s);
    
    fprintf('For season %s:\n', season_names{s});
    
    % fit distributions
    pd_birnbaumsaunders = fitdist(season_bikes, 'BirnbaumSaunders');
    pd_exponential = fitdist(season_bikes, 'Exponential');
    pd_extremevalue = fitdist(season_bikes, 'Extreme Value');
    pd_gamma = fitdist(season_bikes, 'Gamma');
    pd_genextremevalue = fitdist(season_bikes, 'Generalized Extreme Value');
    pd_genpareto = fitdist(season_bikes, 'Generalized Pareto');
    pd_halfnormal = fitdist(season_bikes, 'Half Normal');
    pd_inversegaussian = fitdist(season_bikes, 'Inverse Gaussian');
    pd_logistic = fitdist(season_bikes, 'Logistic');
    pd_loglogistic = fitdist(season_bikes, 'Loglogistic');
    pd_lognormal = fitdist(season_bikes, 'Lognormal');
    pd_nakagami = fitdist(season_bikes, 'Nakagami');
    pd_normal = fitdist(season_bikes, 'Normal');
    pd_rayleigh = fitdist(season_bikes, 'Rayleigh');
    pd_rician = fitdist(season_bikes, 'Rician');
    pd_stable = fitdist(season_bikes, 'Stable');
    pd_tlocationscale = fitdist(season_bikes, 'tlocationscale');
    pd_weibull = fitdist(season_bikes, 'Weibull');
    
    
    % evaluate goodness-of-fit
    
    % The evaluation of goodness of fit was originally going to be conducted
    % using the chi2gof function. However, the p-values were extremely low
    % and no distribution was considered a good enough fit for all the
    % seasons. That is the reason the Kolmogorov-Smirnov test is employed,
    % even though the data is binned.
    % It should be noted that the version this code is developed is the
    % R2016a, so it is possible that the results are affected.
    
    [h_birnbaumsaunders, p_birnbaumsaunders] = kstest(season_bikes, 'CDF', pd_birnbaumsaunders);
    [h_exponential, p_exponential] = kstest(season_bikes, 'CDF', pd_exponential);
    [h_extremevalue, p_extremevalue] = kstest(season_bikes, 'CDF', pd_extremevalue);
    [h_gamma, p_gamma] = kstest(season_bikes, 'CDF', pd_gamma);
    [h_genextremevalue, p_genextremevalue] = kstest(season_bikes, 'CDF', pd_genextremevalue);
    [h_genpareto, p_genpareto] = kstest(season_bikes, 'CDF', pd_genpareto);
    [h_halfnormal, p_halfnormal] = kstest(season_bikes, 'CDF', pd_halfnormal);
    [h_inversegaussian, p_inversegaussian] = kstest(season_bikes, 'CDF', pd_inversegaussian);
    [h_logistic, p_logistic] = kstest(season_bikes, 'CDF', pd_logistic);
    [h_loglogistic, p_loglogistic] = kstest(season_bikes, 'CDF', pd_loglogistic);
    [h_lognormal, p_lognormal] = kstest(season_bikes, 'CDF', pd_lognormal);
    [h_nakagami, p_nakagami] = kstest(season_bikes, 'CDF', pd_nakagami);
    [h_normal, p_normal] = kstest(season_bikes, 'CDF', pd_normal);
    [h_rayleigh, p_rayleigh] = kstest(season_bikes, 'CDF', pd_rayleigh);
    [h_rician, p_rician] = kstest(season_bikes, 'CDF', pd_rician);
    [h_stable, p_stable] = kstest(season_bikes, 'CDF', pd_stable);
    [h_tlocationscale, p_tlocationscale] = kstest(season_bikes, 'CDF', pd_tlocationscale);
    [h_weibull, p_weibull] = kstest(season_bikes, 'CDF', pd_weibull);
    
    % print results
    fprintf('Birnbaum-Saunders distribution: h = %d with p-value = %f\n', h_birnbaumsaunders, p_birnbaumsaunders);
    fprintf('Exponential distribution: h = %d with p-value = %f\n', h_exponential, p_exponential);
    fprintf('Extreme Value distribution: h = %d with p-value = %f\n', h_extremevalue, p_extremevalue);
    fprintf('Gamma distribution: h = %d with p-value = %f\n', h_gamma, p_gamma);
    fprintf('Generalized Extreme Value distribution: h = %d with p-value = %f\n', h_genextremevalue, p_genextremevalue);
    fprintf('Generalized Pareto distribution: h = %d with p-value = %f\n', h_genpareto, p_genpareto);
    fprintf('Half Normal distribution: h = %d with p-value = %f\n', h_halfnormal, p_halfnormal);
    fprintf('Inverse Gaussian distribution: h = %d with p-value = %f\n', h_inversegaussian, p_inversegaussian);
    fprintf('Logistic distribution: h = %d with p-value = %f\n', h_logistic, p_logistic);
    fprintf('Log-logistic distribution: h = %d with p-value = %f\n', h_loglogistic, p_loglogistic);
    fprintf('Log-normal distribution: h = %d with p-value = %f\n', h_lognormal, p_lognormal);
    fprintf('Nakagami distribution: h = %d with p-value = %f\n', h_nakagami, p_nakagami);
    fprintf('Normal distribution: h = %d with p-value = %f\n', h_normal, p_normal);
    fprintf('Rayleigh distribution: h = %d with p-value = %f\n', h_rayleigh, p_rayleigh);
    fprintf('Rician distribution: h = %d with p-value = %f\n', h_rician, p_rician);
    fprintf('Stable distribution: h = %d with p-value = %f\n', h_stable, p_stable);
    fprintf('T-Location-Scale distribution: h = %d with p-value = %f\n', h_tlocationscale, p_tlocationscale);
    fprintf('Weibull distribution: h = %d with p-value = %f\n', h_weibull, p_weibull);
    
    % find best fitting distribution
    p_val = [p_birnbaumsaunders, p_exponential, p_extremevalue, ...
         p_gamma, p_genextremevalue, p_genpareto, p_halfnormal, p_inversegaussian, ...
         p_logistic, p_loglogistic, p_lognormal, p_nakagami, ...
         p_normal, p_rayleigh, p_rician, p_stable, p_tlocationscale, p_weibull];

    h_val = [h_birnbaumsaunders, h_exponential, h_extremevalue, ...
         h_gamma, h_genextremevalue, h_genpareto, h_halfnormal, h_inversegaussian, ...
         h_logistic, h_loglogistic, h_lognormal, h_nakagami, ...
         h_normal, h_rayleigh, h_rician, h_stable, h_tlocationscale, h_weibull];

    distributions = {'Birnbaum-Saunders', 'Exponential', 'Extreme Value', ...
                 'Gamma', 'Generalized Extreme Value', 'Generalized Pareto', 'Half Normal', 'Inverse Gaussian', ...
                 'Logistic', 'Loglogistic', 'Lognormal', 'Nakagami', ...
                 'Normal', 'Rayleigh', 'Rician', 'Stable', 'T-Location-Scale', 'Weibull'};
    
    % find max p-value
    [max_p, max_i] = max(p_val);
    
    % print results of best fit
    if h_val(max_i) == 0
        fprintf('For the season of %s the best fit is the %s distribution.\n', season_names{s}, distributions{max_i});
    else
        fprintf('For the season of %s no distribution is a good enough fit.\n', season_names{s});
        fprintf('But the distribution with highest p-value is %s.\n', distributions{max_i})
    end
    fprintf('\n');
    
    % We plot the distributions with p-values greater than 10^(-5) to better
    % visualize the results. We include p-values that don't meet the qualifying
    % threshold because, as we will see, there is one season without a
    % well-fitting distribution. In this way, we can visualize the data
    % histogram alongside distributions that nearly explained the data but
    % ultimately did not meet the criteria.
    
    % create subplots
    subplot(2, 2, s);
    hold on;
    
    % plot histogram
    histogram(season_bikes, 'Normalization', 'pdf', 'DisplayName', 'Data');
    
    x = linspace(min(season_bikes), max(season_bikes), 100);
    
    % Fit and plot Gamma Distribution
    if p_gamma > 0.00001
        gammaPDF = gampdf(x, pd_gamma.a, pd_gamma.b);
        plot(x, gammaPDF, 'b-', 'LineWidth', 2, 'DisplayName', 'Gamma Fit');
    end
    
    % Fit and plot Exponential Distribution
    if p_exponential > 0.00001
        exponentialPDF = exppdf(x, pd_exponential.mu);
        plot(x, exponentialPDF, 'b--', 'LineWidth', 2, 'DisplayName', 'Exponential Fit');
    end
    
    % Fit and plot Log-Normal Distribution
    if p_lognormal > 0.00001
        lognormalPDF = lognpdf(x, pd_lognormal.mu, pd_lognormal.sigma);
        plot(x, lognormalPDF, 'b:', 'LineWidth', 2, 'DisplayName', 'Log-Normal Fit');
    end
    
    % Fit and plot Generalized Extreme Value Distribution
    if p_genextremevalue > 0.00001
        genextremevaluePDF = pdf(pd_genextremevalue, x);
        plot(x, genextremevaluePDF, 'r-', 'LineWidth', 2, 'DisplayName', 'Generalized Extreme Value Fit');
    end
    
    % Fit and plot Generalized Pareto Distribution
    if p_genpareto > 0.00001
        genparetoPDF = pdf(pd_genpareto, x);
        plot(x, genparetoPDF, 'r--', 'LineWidth', 2, 'DisplayName', 'Generalized Pareto Fit');
    end
    
    % Fit and plot Half-Normal Distribution
    if p_halfnormal > 0.00001
        halfnormalPDF = pdf(pd_halfnormal, x);
        plot(x, halfnormalPDF, 'r:', 'LineWidth', 2, 'DisplayName', 'Half-Normal Fit');
    end
    
    % Fit and plot Logistic Distribution
    if p_logistic > 0.00001
        logisticPDF = pdf(pd_logistic, x);
        plot(x, logisticPDF, 'g:', 'LineWidth', 2, 'DisplayName', 'Logistic Fit');
    end
    
    % Fit and plot Log-Logistic Distribution
    if p_loglogistic > 0.00001
        loglogisticPDF = pdf(pd_loglogistic, x);
        plot(x, loglogisticPDF, 'g-', 'LineWidth', 2, 'DisplayName', 'Log-Logistic Fit');
    end
    
    % Fit and plot Nakagami Distribution
    if p_nakagami > 0.00001
        nakagamiPDF = pdf(pd_nakagami, x);
        plot(x, nakagamiPDF, 'c-', 'LineWidth', 2, 'DisplayName', 'Nakagami Fit');
    end
    
    % Fit and plot Rayleigh Distribution
    if p_rayleigh > 0.00001
        rayleighPDF = pdf(pd_rayleigh, x);
        plot(x, rayleighPDF, 'c--', 'LineWidth', 2, 'DisplayName', 'Rayleigh Fit');
    end
    
    % Fit and plot Rician Distribution
    if p_rician > 0.00001
        ricianPDF = pdf(pd_rician, x);
        plot(x, ricianPDF, 'c:', 'LineWidth', 2, 'DisplayName', 'Rician Fit');
    end
    
    % Fit and plot Stable Distribution
    if p_stable > 0.00001
        stablePDF = pdf(pd_stable, x);
        plot(x, stablePDF, 'g--', 'LineWidth', 2, 'DisplayName', 'Stable Fit');
    end
    
    % Fit and plot T-Location-Scale Distribution
    if p_tlocationscale > 0.00001
        tlocationscalePDF = pdf(pd_tlocationscale, x);
        plot(x, tlocationscalePDF, 'm--', 'LineWidth', 2, 'DisplayName', 'T-Location-Scale Fit');
    end
    
    % Fit and plot Weibull Distribution
    if p_weibull > 0.00001
        weibullPDF = pdf(pd_weibull, x);
        plot(x, weibullPDF, 'm-', 'LineWidth', 2, 'DisplayName', 'Weibull Fit');
    end
    
    xlabel('Number of Bikes Rented');
    ylabel('Probalility Density');
    title(['Season ', season_names{s}]);
    legend('show');
end

hold off;

% Upon running the program, we observe that all p-values are relatively low.
% However, there are some distributions that we cannot reject as fitting
% the histogram data.

% For the Winter season, the only p-value high enough to not reject the
% hypothesis was for the Weibull distribution, with a p-value of 0.138541.
% While several distributions appear close to fitting the data, as seen from
% the plots, only the Weibull distribution has a sufficiently high p-value.

% For the Spring season, only the Nakagami distribution has a sufficiently
% high p-value of 0.088690. As seen in the plot, the Spring season has the
% fewest drawn distributions, indicating that most distributions failed to
% surpass even the 10^(-5) threshold we established.

% The Summer season is similar to Spring, as once again, only the Nakagami
% distribution cannot be rejected, with a p-value of 0.064111. As we can
% see from the plots, , a few more distributions come close to describing
% the data, but ultimately fall short.

% The Autumn season is peculiar, as no distribution comes close to having
% an acceptable p-value. The closest is the Nakagami distribution, with a
% p-value of 0.003115, which is still one order of magnitude lower than the
% desired threshold. This failure to identify a distribution that accurately
% explains the data is highlighted in the plot, where we can see that even
% the highest p-value distributions do not align with the histogram data.
