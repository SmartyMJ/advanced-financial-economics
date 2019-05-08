
function output = descriptiveStats(input)

%This provides the Excel descriptive stats for a given vector
    % Stats are given in the following order: 
    % mean, std error, median, mode, std. deviation, sample variance, kurtosis, skewness, range, min, max, sum, count

stderror = std(input) / sqrt(length(input));
output = [mean(input); stderror; median(input); mode(input); std(input); var(input); kurtosis(input); skewness(input); range(input); min(input); max(input); sum(input); length(input)];