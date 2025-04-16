% Analush kai Epeksergasia Dedomenwn
% Zhthma 2 
% Papadopoulou Martha

clear; clc; close all;

% import data
data = importdata('SeoulBike.xlsx');
bikes = data.data(:,1); %rented bike count
season = data.data(:,end-1);
season_names = {'Winter', 'Spring', 'Summer', 'Autumn'};
colors = {'cyan', 'magenta', 'yellow', 'red'};

figure(1); % create a figure
iter_count = 1; % count of iteration
bin_count = 20; % number of bins
percentage = zeros(6,1); % initialize vector for percentage values

fprintf('The percentage of times the first season''s data explains the second season''s data is:\n');

for i = 1:4
    for j = (i+1):4
        % find indices for correct seasons
        index_1 = (season == i);
        index_2 = (season == j);
        
        % access bike data for correct seasons
        season_bikes_1 = bikes(index_1);
        season_bikes_2 = bikes(index_2);
                
        % calculate bin edges
        % The reason for manually calculating the bin edges is to ensure
        % that all histograms use the same bin edges, which is necessary
        % for accurately calculating the chi-square goodness of fit.
        bin_edges = linspace(min([season_bikes_1; season_bikes_2]), max([season_bikes_1; season_bikes_2]), bin_count + 1);
        
        % plot histograms for both seasons
        subplot(2, 3, iter_count);
        histogram(season_bikes_1, 'BinEdges', bin_edges, 'FaceColor', colors{i}, 'FaceAlpha', 0.5);
        hold on;
        histogram(season_bikes_2, 'BinEdges', bin_edges, 'FaceColor', colors{j}, 'FaceAlpha', 0.5);
        
        xlabel('Rented Bike Count');
        ylabel('Frequency');
        title([season_names{i}, ' vs ', season_names{j}]);
        legend(season_names{i}, season_names{j});
        
        n = 100; % number of random observetions
        M = 100; % number of repetitions
        count = 0; % count of times the null hypothesis is not rejected
        
        for k = 1:M
            % random indices for each season
            rand_index_1 = randperm(length(season_bikes_1), n);
            rand_index_2 = randperm(length(season_bikes_2), n);
            
            % random bike data for each season
            rand_bikes_1 = season_bikes_1(rand_index_1);
            rand_bikes_2 = season_bikes_2(rand_index_2);
            
            % expected counts from season_bikes_1
            expected_counts = histcounts(rand_bikes_1, bin_edges);
            % observed counts from season_bikes_2
            observed_counts = histcounts(rand_bikes_2, bin_edges);
            
            % filter zero data to avoid divition with 0
            non_zero_indices = expected_counts > 0;
            observed_counts = observed_counts(non_zero_indices);
            expected_counts = expected_counts(non_zero_indices);
            
            % scale expected counts
            % This extra step is needed because the counts of each season
            % are different. In order to compare the seasons, we need to
            % calculate the percentage of counts in the expected data and
            % then convert this to match the total counts of the observed
            % data.
            expected_prob = expected_counts / sum(expected_counts);
            expected_counts_scaled = expected_prob * sum(observed_counts);
            
            % manually calculate the chi-square statistic and p-value
            chi2_stat = sum(((observed_counts - expected_counts_scaled).^2) ./ expected_counts_scaled);
            df = length(expected_counts_scaled) - 1; % degrees of freedom
            p_value = 1 - chi2cdf(chi2_stat, df);
            
            % check if the null hypothesis is not rejected
            if p_value > 0.05
                count = count + 1;
            end
            
        end
        
        % calculate percentage of times the null hypothesis was not rejected
        percentage(iter_count) = count/M;
        
        % print results
        fprintf('For the seasons %s and %s: %.2f%%\n', season_names{i}, season_names{j}, percentage(iter_count)*100);
        
        iter_count = iter_count + 1; % for next iteration
        
    end
    
end

hold off;

% From the plots, it is evident that the distributions for the seasons do
% not closely follow each other, especially with Winter being the most
% distinct. The plots show that the Winter season differs significantly from
% the others, with the majority of its counts concentrated in the first
% bins. This observation us further supported by the calculation results.

% Since the results depend on the random selection of data, they vary
% between runs. To gain a broader understanding of the potential outcomes,
% we run the program multiple times.

% For Winter as the expected distribution and Spring as the observed data,
% the percentage at which the observed data aligns with the expected
% distribution is between 0% and 5%. 

% For Winter as the first season and Summer as the second, the percentage
% of alignment is between 0% and 1%, with 0% occuring most frequently.

% In a similar manner, Winter's distribution explains Autumn's data by only
% 0% to 2%.

% For the remaining seasons' comparisons, the percentages are slightly
% higher but still not substantial.

% Spring's distribution aligns with Summer's data by a percentage ranging
% from 0% to 6%.

% Similarly, when comparing Spring as the first season and Autumn as the 
% second, the percentage alignment ranges from 5% to 18%.

% Lastly, the alignment between Summer's distribution and Autumn's data is
% relatively high, ranging from 16% to 26%.

% It is important to note that these results are closely tied to the number
% of bins used to categorize the data. Any changes in the binning would
% directly impact the outcomes.
