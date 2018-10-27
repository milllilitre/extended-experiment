function DrawFrequencyComparison()

% ava_dists_of_methods = rand(1, 16) * 0.07 + 0.545;
% variance_of_methods = rand(1, 16) * 0.07 + 0.08;
load('frequencyComparison.mat');
figure(); hold on;
    bar_obj = bar(ava_dists_of_methods);
    bar_obj.FaceColor = 'flat';
    bar_obj.CData(:,:) = repmat([.5 0 .5], 16, 1);
    errorbar(ava_dists_of_methods, variance_of_methods, 'k', 'lineStyle', 'none');
    
ylabel('Error(m)');
xlabel('Channel');
xticks(1:16);
end

