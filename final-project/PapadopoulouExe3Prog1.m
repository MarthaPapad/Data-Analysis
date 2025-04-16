% Analush kai Epeksergasia Dedomenwn
% Zhthma 3
% Papadopoulou Martha
% AEM: 4438

clear; clc; close all;

% import data
data = importdata('SeoulBike.xlsx');
bikes = data.data(:,1); %rented bike count
hours = data.data(:,2);
season = data.data(:,end-1);
season_names = {'Winter', 'Spring', 'Summer', 'Autumn'};

% filter data for specific season
spec_season = 2; % spring season (can be changed)
season_bikes = bikes(season == spec_season);
season_hours = hours(season == spec_season);

% find the indices where the days start
days_index_logical = (season_hours == 0);
days_index = find(days_index_logical);

% initialize a matrix for the hours' difference
days_num = length(days_index);
difference = zeros(days_num, 276);

% initialize matrices for mean differences and t-test results
mean_hour_diff = NaN(24);
h_hour = NaN(24);

% counter of unique hour pairs
counter = 1;

% calculate difference for all days of the season
for i = 1:24
    for j = i+1:24
        for d =1:days_num
            % select specific day data
            spec_day_index = days_index(d);
            spec_day_bikes = season_bikes(spec_day_index:spec_day_index+23);
            
            % calculate differences
            difference(d, counter) = spec_day_bikes(i) - spec_day_bikes(j);
        end
        
        % calculate mean difference and t-test results of an hour pair
        mean_hour_diff(i,j) = mean(difference(:, counter));
        [h_hour(i,j), ~] = ttest(difference(:,counter), 0);
        
        counter = counter + 1;
    end
end


% visualize the results using imagesc (heatmap is unavailable in the 2016a version)
% imagesc for mean differences of each hour pair
figure(1);
imagesc(mean_hour_diff);
colormap(spring); 
colorbar;
title(['Mean Differences of Bike Count Between Hours for ' season_names{spec_season}, ' Season']);
xlabel('Hour');
ylabel('Hour');
set(gca, 'XTick', 1:24, 'YTick', 1:24);

% add values on the squares
for i = 1:24
    for j = 1:24
        text(j, i, num2str(mean_hour_diff(i, j), '%.0f'), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'black');
    end
end

% imagesc for the t-test results
figure(2);
imagesc(h_hour);
colormap([1 1 0.5; 1 0.5 1]); % yellow for 0, magenta for 1
colorbar;
title(['T-test Results for the Null Hypothesis Between Hours for ' season_names{spec_season}, ' Season']);
xlabel('Hour');
ylabel('Hour');
set(gca, 'XTick', 1:24, 'YTick', 1:24);

% add values on the squares
for i = 1:24
    for j = 1:24
        text(j, i, num2str(h_hour(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'black');
    end
end

% The results discussed come from the spring season. Although the mean
% differences heatmap remains fairly consistent across all seasons, the
% t-test results show more variation. It is important to note that only the
% differences for the unique 276 pairs of hours were calculated.  As a
% result, the heatmap from the diagonal and downward contains NaN values, since
% these represent the same comparisons as those above the diagonal, but
% with opposite signs, making further calculation unnecessary.

% From the heatmap of mean differences in bike counts per hour pair, we
% observe significant and frequent deviations from 0. The largest
% differences occur between the 19th hour of the day and the first 12 hours,
% as well as the 24th hour. This heatmap further emphasizes the fact that rented
% bike counts vary with the hour of the day. It is common for some hours to
% have higher bike rentals than others, so the mean difference in bike
% counts is rarely 0.

% The second heatmap shows the number of hour pairs for which we cannot 
% reject the hypothesis that their mean difference is close to 0. As
% previously mentioned, the NaN values along the diagonal and below are
% excluded from the analysis. Ignoring these, we observe that only a small
% number of hour pairs have h = 0. Out of the 276 pairs, only 21 cannot be
% rejected.
