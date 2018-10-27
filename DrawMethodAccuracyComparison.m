function DrawMethodAccuracyComparison(avaDists, variance, medDists)

variance_of_methods = zeros(1, 6);
ava_dists_of_methods = [1.504 0.90 0.7 0.62 0.687 0.55];
figure();
bar_obj = bar(ava_dists_of_methods, 0.5);
bar_obj.FaceColor = 'flat';
ylabel('Average Error(m)');
xlabel('Methods');
xticks(1:6);
xticklabels({'LANDMARC','Li''s method', 'TagTrack','MaTrack','IMPSO-FNN', 'WallSense'});
bar_obj.CData(:,:) = repmat([.5 0 .5], 6, 1);


end