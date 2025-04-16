% Analush kai Epeksergasia Dedomenwn
% Zhthma 4
% Papadopoulou Martha

clear; clc; close all;

% import data
data = importdata('SeoulBike.xlsx');
bikes = data.data(:,1); %rented bike count
hours = data.data(:,2);
season = data.data(:,end-1);
holiday = data.data(:,end);
season_names = {'Winter', 'Spring', 'Summer', 'Autumn'};

% divide data whether it is holiday or not
holiday_bikes = bikes(holiday == 1);
non_holiday_bikes = bikes(holiday == 0);

% calculate the medians for all seasons included
mu_1 = median(holiday_bikes);
mu_0 = median(non_holiday_bikes);

% select specific season
spec_season = 3; % summer season (can be changed)

% initialize median and confidence interval matrices
M = 1000; % number of iterations for bootstrap
median_mu_1 = zeros(24, 1);
median_mu_0 = zeros(24, 1);
ci_mu_1 = zeros(24, 2);
ci_mu_0 = zeros(24, 2);
a = 0.05; % 95% condfidence level

for h = 1:24
    % select specific hour
    spec_hour = h-1; 
    
    % filter data for those specific values 
    spec_hol_bikes = bikes(holiday == 1 & season == spec_season & hours == spec_hour);
    spec_non_hol_bikes = bikes(holiday == 0 & season == spec_season & hours == spec_hour);
    
    % bootstrap
    n_1 = length(spec_hol_bikes);
    n_0 = length(spec_non_hol_bikes);
    boot_mu_1 = zeros(M,1);
    boot_mu_0 = zeros(M,1);
    
    for i = 1:M
        % generate random indices
        index_1 = unidrnd(n_1, n_1, 1);
        index_0 = unidrnd(n_0, n_0, 1);
        
        % access random indices data
        rand_hol_bikes = spec_hol_bikes(index_1);
        rand_non_hol_bikes = spec_non_hol_bikes(index_0);
        
        % calculate the medians
        boot_mu_1(i) = median(rand_hol_bikes);
        boot_mu_0(i) = median(rand_non_hol_bikes);
    end
    
    % calculate medians of medians
    median_mu_1(h) = median(boot_mu_1);
    median_mu_0(h) = median(boot_mu_0);
    
    % calculate confidence interval
    sorted_mu_1 = sort(boot_mu_1);
    sorted_mu_0 = sort(boot_mu_0);
    
    % find lower and upper index
    lower_index = ceil(M * a/2);
    upper_index = floor(M * (1 - a/2));
    
    % add confidence intervals' values in their matrices
    ci_mu_1(h, :) = [sorted_mu_1(lower_index), sorted_mu_1(upper_index)];
    ci_mu_0(h, :) = [sorted_mu_0(lower_index), sorted_mu_0(upper_index)];
    
end

% plot for holidays
figure(1);
errorbar(0:23, median_mu_1, median_mu_1 - ci_mu_1(:,1), ci_mu_1(:,2) - median_mu_1, 'o-', 'LineWidth', 1.5);
hold on;
plot(-1:24, mu_1 * ones(1, 26), 'r--', 'LineWidth', 1.5); 
xlabel('Hour of the day');
ylabel('Median Bike Rentals (Holiday)');
title(['95% Bootstrap CI for Median Bike Rentals (Holiday) for ', season_names{spec_season}]);
xlim([-1 24]);
legend('Hourly Medians with 95% CI', 'Overall Median \mu_1', 'Location', 'northwest');
grid on;

% plot for non-holidays
figure(2);
errorbar(0:23, median_mu_0, median_mu_0 - ci_mu_0(:,1), ci_mu_0(:,2) - median_mu_0, 'o-', 'LineWidth', 1.5);
hold on;
plot(-1:24, mu_0 * ones(1, 26), 'r--', 'LineWidth', 1.5);
xlabel('Hour of the day');
ylabel('Median Bike Rentals (Non-Holiday)');
title(['95% Bootstrap CI for Median Bike Rentals (Non-Holiday) for ', season_names{spec_season}]);
xlim([-1 24]);
legend('Hourly Medians with 95% CI', 'Overall Median \mu_0', 'Location', 'northwest');
grid on;

hold off;

% The results discussed here are for the summer season, though all seasons
% show similar outcomes. The season under study can be changed, and the
% respective results will be displayed in the plots. Across all seasons,
% the number of holiday days is very low, with summer having the fewest, 2
% days. This is the reason why the confidence intervals for holidays are
% the widest, since the bootstrap is sampling from just 2 values.

% For the holiday plot, we observe that during the first 7 hours of the day, 
% the median bike rentals have relatively narrow confidence intervals.
% However, for nearly all other hours, the confidence intervals are
% extremely wide. This is expected because, as mentioned earlier, the
% bootstrap process only has two numbers to sample from. The overall median
% mu_1 is quite low at 259, with only the 4th, 5th, and 6th hours coming
% close to this value. This makes sense, as the overall median is calculated
% from all days, hours, and seasons, so it naturally differs from individual
% hourly medians due to the high variability in bike rentals, which are
% influenced by these factors. The hourly medians continuously drop at first
% until the 5th hour, then begin to rise, peaking at the 20th hour, before
% they start to fall again.

% On the other hand, non-holiday days have a much larger sample size of 90
% days for the summer season. So, as expected, the confidence intervals are
% significantly narrower than those for holidays. Each hour of the day has
% a distinct  median rental value, except for the hours between 10th and
% 14th hour of the day, where the median remains relatively stable. The
% hourly medians initially drop, then rise and peak at the 8th hour, before
% dropping again and remaining steady from the 10th to the 14th hour.
% They rise once more, peaking at the 18th hour, and then drop again. The
% overall median for non-holiday days mu_0 is moderately low at 561.
% However, unlike in the holiday case, there are hours at which the median
% rentals fall below the overall value, between the 10th and 15th hour. The
% only hourly medians that are close to the overall median are the 2nd and
% 6th hour. 
