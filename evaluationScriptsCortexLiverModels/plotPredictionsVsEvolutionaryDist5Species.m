% Get the distances of glires from Mus musculus
distanceFromMouseGlires = [28.6
55
71
73
70
73
73
73
73
33
33
73
73
73
73
73
70
70
73
33
33
73
71
71
73
73
73
71
55
82
71
28.6
33
33
7.04
0
8.3
3.07
71
73
45
82
73
33
33
82
70
33
73
28.6
20.9
33
71
73
71
55
];

% Plot evolutionary distance vs. predictions for non-enhancer orthologs of
% enhancers negative set for cortex mouse-only model
predictedProba = importdata('Pfenning_bulk_Ctx_nonCDS_enhancerShort_noRepeats_test_andNeg_conv5MuchMuchMoreFiltNoInterPoolTwoDenseLargeDenseVeryHighMomLowDropL2SmallBatchPretrainBal_predictedProba.txt');
mean(predictedProba(15925:26296)) % 0.2506
mousePeakPredictionsTestGliresCortex = importdata('Pfenning_bulk_Ctx_nonCDS_enhancerShort_test_peakPredictionsMouseOnlyGlires_noRepeats.txt');
mousePeakPredictionsTestGliresCortex.data(find(mousePeakPredictionsTestGliresCortex.data == -1)) = NaN;
meanMousePeakPredictionsTestGliresCortex = nanmean(mousePeakPredictionsTestGliresCortex.data)';
mousePeakPredictionsIndicatorTestGliresCortex = zeros(size(mousePeakPredictionsTestGliresCortex.data, 1), 1);
for i=1:size(mousePeakPredictionsTestGliresCortex.data, 1)
if length(find(isnan(mousePeakPredictionsTestGliresCortex.data(i, :)))) < 14
mousePeakPredictionsIndicatorTestGliresCortex(i) = 1;
end
end
sum(mousePeakPredictionsIndicatorTestGliresCortex)
mousePeakPredictionsFiltTestGliresCortex = mousePeakPredictionsTestGliresCortex.data(find(mousePeakPredictionsIndicatorTestGliresCortex==1), :);
meanMousePeakPredictionsFiltTestGliresCortex = nanmean(mousePeakPredictionsFiltTestGliresCortex)';
[c, p] = corr(meanMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires) % c = -0.9308, p = 2.8757e-25
[c, p] = corr(meanMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = -0.6308, p = 1.8729e-07
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMouseCortex=sortrows([distanceFromMouseGlires meanMousePeakPredictionsFiltTestGliresCortex]);
sortedMouseCortex_x = sortedMouseCortex(:,1);
sortedMouseCortex_y = sortedMouseCortex(:,2);
yy1MouseCortex = smooth(sortedMouseCortex_x, sortedMouseCortex_y, 0.9, 'loess');
plot(sortedMouseCortex_x, sortedMouseCortex_y, 'k.', sortedMouseCortex_x, yy1MouseCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.2506, 0.2506], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
stdMousePeakPredictionsFiltTestGliresCortex = nanstd(mousePeakPredictionsFiltTestGliresCortex)';
[c, p] = corr(stdMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires) % c = 0.7634, p = 7.8403e-12
[c, p] = corr(stdMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = 0.7586, p = 1.2628e-11
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMouseCortex=sortrows([distanceFromMouseGlires stdMousePeakPredictionsFiltTestGliresCortex]);
sortedMouseCortexStd_x = sortedMouseCortex(:,1);
sortedMouseCortexStd_y = sortedMouseCortex(:,2);
yy1MouseStdCortex = smooth(sortedMouseCortexStd_x, sortedMouseCortexStd_y, 0.9, 'loess');
plot(sortedMouseCortexStd_x, sortedMouseCortexStd_y, 'k.', sortedMouseCortexStd_x, yy1MouseStdCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for non-enhancer orthologs of
% enhancers negative set for cortex mouse-only model with N filtering
mousePeakPredictionsTestGliresCortex = importdata('Pfenning_bulk_Ctx_nonCDS_enhancerShort_test_peakPredictionsMouseOnlyMaxNs0.05Glires_noRepeats.txt');
mousePeakPredictionsTestGliresCortex.data(find(mousePeakPredictionsTestGliresCortex.data == -1)) = NaN;
mousePeakPredictionsIndicatorTestGliresCortex = zeros(size(mousePeakPredictionsTestGliresCortex.data, 1), 1);
for i=1:size(mousePeakPredictionsTestGliresCortex.data, 1)
if length(find(isnan(mousePeakPredictionsTestGliresCortex.data(i, :)))) < 14
mousePeakPredictionsIndicatorTestGliresCortex(i) = 1;
end
end
sum(mousePeakPredictionsIndicatorTestGliresCortex)
mousePeakPredictionsFiltTestGliresCortex = mousePeakPredictionsTestGliresCortex.data(find(mousePeakPredictionsIndicatorTestGliresCortex==1), :);
meanMousePeakPredictionsFiltTestGliresCortex = nanmean(mousePeakPredictionsFiltTestGliresCortex)';
[c, p] = corr(meanMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires) % c = -0.9334, p = 1.0664e-25
[c, p] = corr(meanMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = -0.6168, p = 4.1494e-07
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMouseCortex=sortrows([distanceFromMouseGlires meanMousePeakPredictionsFiltTestGliresCortex]);
sortedMouseCortex_x = sortedMouseCortex(:,1);
sortedMouseCortex_y = sortedMouseCortex(:,2);
yy1MouseCortex = smooth(sortedMouseCortex_x, sortedMouseCortex_y, 0.9, 'loess');
plot(sortedMouseCortex_x, sortedMouseCortex_y, 'k.', sortedMouseCortex_x, yy1MouseCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.2506, 0.2506], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
stdMousePeakPredictionsFiltTestGliresCortex = nanstd(mousePeakPredictionsFiltTestGliresCortex)';
[c, p] = corr(stdMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires) % c = 0.7712, p = 3.5351e-12
[c, p] = corr(stdMousePeakPredictionsFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = 0.7764, p = 2.0432e-12
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMouseCortex=sortrows([distanceFromMouseGlires stdMousePeakPredictionsFiltTestGliresCortex]);
sortedMouseCortexStd_x = sortedMouseCortex(:,1);
sortedMouseCortexStd_y = sortedMouseCortex(:,2);
yy1MouseStdCortex = smooth(sortedMouseCortexStd_x, sortedMouseCortexStd_y, 0.9, 'loess');
plot(sortedMouseCortexStd_x, sortedMouseCortexStd_y, 'k.', sortedMouseCortexStd_x, yy1MouseStdCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for non-enhancer orthologs of
% enhancers negative set for cortex multi-species model
predictedProbaMulti = importdata('Pfenning_bulk_Ctx_nonCDS_enhancerShort_test_summitPlusMinus250bp_andMacaque_andRat_andBat_andNeg_predictedProba.txt');
mean(predictedProbaMulti(41641:79278)) % 0.2599
mousePeakPredictionsMultiTestGliresCortex = importdata('Pfenning_bulk_Ctx_nonCDS_enhancerShort_test_peakPredictionsMultiSpeciesGlires_noRepeats.txt');
mousePeakPredictionsMultiTestGliresCortex.data(find(mousePeakPredictionsMultiTestGliresCortex.data == -1)) = NaN;
mousePeakPredictionsMultiIndicatorTestGliresCortex = zeros(size(mousePeakPredictionsMultiTestGliresCortex.data, 1), 1);
for i=1:size(mousePeakPredictionsMultiTestGliresCortex.data, 1)
if length(find(isnan(mousePeakPredictionsMultiTestGliresCortex.data(i, :)))) < 14
mousePeakPredictionsMultiIndicatorTestGliresCortex(i) = 1;
end
end
sum(mousePeakPredictionsMultiIndicatorTestGliresCortex)
mousePeakPredictionsMultiFiltTestGliresCortex = mousePeakPredictionsMultiTestGliresCortex.data(find(mousePeakPredictionsMultiIndicatorTestGliresCortex==1), :);
meanMousePeakPredictionsMultiFiltTestGliresCortex = nanmean(mousePeakPredictionsMultiFiltTestGliresCortex)';
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires) % c = -0.9676, p = 6.0776e-34
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = -0.7897, p = 4.6882e-13
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiCortex=sortrows([distanceFromMouseGlires meanMousePeakPredictionsMultiFiltTestGliresCortex]);
sortedMultiCortex_x = sortedMultiCortex(:,1);
sortedMultiCortex_y = sortedMultiCortex(:,2);
yy1MultiCortex = smooth(sortedMultiCortex_x, sortedMultiCortex_y, 0.9, 'loess');
plot(sortedMultiCortex_x, sortedMultiCortex_y, 'k.', sortedMultiCortex_x, yy1MultiCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.2599, 0.2599], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
stdMousePeakPredictionsMultiFiltTestGliresCortex = nanstd(mousePeakPredictionsMultiFiltTestGliresCortex)';
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires) % c = 0.8012, p = 1.2092e-13
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = 0.6556, p = 4.1543e-08
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiCortexStd=sortrows([distanceFromMouseGlires stdMousePeakPredictionsMultiFiltTestGliresCortex]);
sortedMultiCortexStd_x = sortedMultiCortexStd(:,1);
sortedMultiCortexStd_y = sortedMultiCortexStd(:,2);
yy1MultiStdCortex = smooth(sortedMultiCortexStd_x, sortedMultiCortexStd_y, 0.9, 'loess');
plot(sortedMultiCortexStd_x, sortedMultiCortexStd_y, 'k.', sortedMultiCortexStd_x, yy1MultiStdCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for non-enhancer orthologs of
% enhancers negative set for cortex multi-species model with N filtering
mousePeakPredictionsMultiTestGliresCortex = importdata('Pfenning_bulk_Ctx_nonCDS_enhancerShort_test_peakPredictionsMultiSpeciesMaxNs0.05Glires_noRepeats.txt');
mousePeakPredictionsMultiTestGliresCortex.data(find(mousePeakPredictionsMultiTestGliresCortex.data == -1)) = NaN;
meanMousePeakPredictionsMultiTestGliresCortex = nanmean(mousePeakPredictionsMultiTestGliresCortex.data)';
mousePeakPredictionsMultiIndicatorTestGliresCortex = zeros(size(mousePeakPredictionsMultiTestGliresCortex.data, 1), 1);
for i=1:size(mousePeakPredictionsMultiTestGliresCortex.data, 1)
if length(find(isnan(mousePeakPredictionsMultiTestGliresCortex.data(i, :)))) < 14
mousePeakPredictionsMultiIndicatorTestGliresCortex(i) = 1;
end
end
sum(mousePeakPredictionsMultiIndicatorTestGliresCortex)
mousePeakPredictionsMultiFiltTestGliresCortex = mousePeakPredictionsMultiTestGliresCortex.data(find(mousePeakPredictionsMultiIndicatorTestGliresCortex==1), :);
meanMousePeakPredictionsMultiFiltTestGliresCortex = nanmean(mousePeakPredictionsMultiFiltTestGliresCortex)';
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires) % c = -0.9680, p = 4.1928e-34
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = -0.7810, p = 1.8729e-07
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiCortex=sortrows([distanceFromMouseGlires meanMousePeakPredictionsMultiFiltTestGliresCortex]);
sortedMultiCortex_x = sortedMultiCortex(:,1);
sortedMultiCortex_y = sortedMultiCortex(:,2);
yy1MultiCortex = smooth(sortedMultiCortex_x, sortedMultiCortex_y, 0.9, 'loess');
plot(sortedMultiCortex_x, sortedMultiCortex_y, 'k.', sortedMultiCortex_x, yy1MultiCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.2599, 0.2599], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
stdMousePeakPredictionsMultiFiltTestGliresCortex = nanstd(mousePeakPredictionsMultiFiltTestGliresCortex)';
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires) % c = 0.8124, p = 2.9356e-14
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresCortex, distanceFromMouseGlires, 'type', 'Spearman') % c = 0.6915, p = 3.6176e-09
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiCortexStd=sortrows([distanceFromMouseGlires stdMousePeakPredictionsMultiFiltTestGliresCortex]);
sortedMultiCortexStd_x = sortedMultiCortexStd(:,1);
sortedMultiCortexStd_y = sortedMultiCortexStd(:,2);
yy1MultiStdCortex = smooth(sortedMultiCortexStd_x, sortedMultiCortexStd_y, 0.9, 'loess');
plot(sortedMultiCortexStd_x, sortedMultiCortexStd_y, 'k.', sortedMultiCortexStd_x, yy1MultiStdCortex, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Get the distances of glires from Mus musculus for alternative species
% ordering
distanceFromMouseGlires = [0
28.6
55
71
73
70
73
73
73
73
33
33
73
73
73
73
73
70
70
73
33
33
73
71
71
73
73
73
71
55
82
71
28.6
33
33
7.04
8.3
3.07
71
73
45
82
73
33
33
82
70
33
73
28.6
20.9
33
71
73
71
55
];

% Plot evolutionary distance vs. predictions for non-enhancer orthologs of
% enhancers negative set for liver multi-species model
predictedProbaMultiLiver = importdata('idr.optimal_peak.inLiuAll_nonCDS_enhancerShort_summitPlusMinus250bp_andMacaque_andRat_andCow_andPig_test_andNeg_predictedProba.txt');
mean(predictedProbaMultiLiver(24535:61250)) % 0.2271
mousePeakPredictionsMultiTestGliresLiver = importdata('idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort_test_peakPredictionsFiveSpeciesGlires.txt');
mousePeakPredictionsMultiTestGliresLiver.data(find(mousePeakPredictionsMultiTestGliresLiver.data == -1)) = NaN;
meanMousePeakPredictionsMultiTestGliresLiver = nanmean(mousePeakPredictionsMultiTestGliresLiver.data)';
mousePeakPredictionsMultiIndicatorTestGliresLiver = zeros(size(mousePeakPredictionsMultiTestGliresLiver.data, 1), 1);
for i=1:size(mousePeakPredictionsMultiTestGliresLiver.data, 1)
if length(find(isnan(mousePeakPredictionsMultiTestGliresLiver.data(i, :)))) < 14
mousePeakPredictionsMultiIndicatorTestGliresLiver(i) = 1;
end
end
sum(mousePeakPredictionsMultiIndicatorTestGliresLiver)
mousePeakPredictionsMultiFiltTestGliresLiver = mousePeakPredictionsMultiTestGliresLiver.data(find(mousePeakPredictionsMultiIndicatorTestGliresLiver==1), :);
meanMousePeakPredictionsMultiFiltTestGliresLiver = nanmean(mousePeakPredictionsMultiFiltTestGliresLiver)';
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires) % c = -0.9224, p = 5.6476e-24
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires, 'type', 'Spearman') % c = -0.6225, p = 3.0172e-07
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiLiver=sortrows([distanceFromMouseGlires meanMousePeakPredictionsMultiFiltTestGliresLiver]);
sortedMultiLiver_x = sortedMultiLiver(:,1);
sortedMultiLiver_y = sortedMultiLiver(:,2);
yy1MultiLiver = smooth(sortedMultiLiver_x, sortedMultiLiver_y, 0.9, 'loess');
plot(sortedMultiLiver_x, sortedMultiLiver_y, 'k.', sortedMultiLiver_x, yy1MultiLiver, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.2271, 0.2271], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
stdMousePeakPredictionsMultiFiltTestGliresLiver = nanstd(mousePeakPredictionsMultiFiltTestGliresLiver)';
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires) % c = 0.3895, p = 0.0030
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires, 'type', 'Spearman') % c = 0.1546, p = 0.2554
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiLiverStd=sortrows([distanceFromMouseGlires stdMousePeakPredictionsMultiFiltTestGliresLiver]);
sortedMultiLiverStd_x = sortedMultiLiverStd(:,1);
sortedMultiLiverStd_y = sortedMultiLiverStd(:,2);
yy1MultiStdLiver = smooth(sortedMultiLiverStd_x, sortedMultiLiverStd_y, 0.9, 'loess');
plot(sortedMultiLiverStd_x, sortedMultiLiverStd_y, 'k.', sortedMultiLiverStd_x, yy1MultiStdLiver, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for non-enhancer orthologs of
% enhancers negative set for liver multi-species model with N filtering
mousePeakPredictionsMultiTestGliresLiver = importdata('idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort_test_peakPredictionsFiveSpeciesMaxNs0.05Glires.txt');
mousePeakPredictionsMultiTestGliresLiver.data(find(mousePeakPredictionsMultiTestGliresLiver.data == -1)) = NaN;
meanMousePeakPredictionsMultiTestGliresLiver = nanmean(mousePeakPredictionsMultiTestGliresLiver.data)';
mousePeakPredictionsMultiIndicatorTestGliresLiver = zeros(size(mousePeakPredictionsMultiTestGliresLiver.data, 1), 1);
for i=1:size(mousePeakPredictionsMultiTestGliresLiver.data, 1)
if length(find(isnan(mousePeakPredictionsMultiTestGliresLiver.data(i, :)))) < 14
mousePeakPredictionsMultiIndicatorTestGliresLiver(i) = 1;
end
end
sum(mousePeakPredictionsMultiIndicatorTestGliresLiver)
mousePeakPredictionsMultiFiltTestGliresLiver = mousePeakPredictionsMultiTestGliresLiver.data(find(mousePeakPredictionsMultiIndicatorTestGliresLiver==1), :);
meanMousePeakPredictionsMultiFiltTestGliresLiver = nanmean(mousePeakPredictionsMultiFiltTestGliresLiver)';
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires) % c = -0.9381, p = 1.5784e-26
[c, p] = corr(meanMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires, 'type', 'Spearman') % c = -0.6146, p = 4.6762e-07
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiLiver=sortrows([distanceFromMouseGlires meanMousePeakPredictionsMultiFiltTestGliresLiver]);
sortedMultiLiver_x = sortedMultiLiver(:,1);
sortedMultiLiver_y = sortedMultiLiver(:,2);
yy1MultiLiver = smooth(sortedMultiLiver_x, sortedMultiLiver_y, 0.9, 'loess');
plot(sortedMultiLiver_x, sortedMultiLiver_y, 'k.', sortedMultiLiver_x, yy1MultiLiver, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.2271, 0.2271], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
stdMousePeakPredictionsMultiFiltTestGliresLiver = nanstd(mousePeakPredictionsMultiFiltTestGliresLiver)';
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires) % c = 0.4666, p = 2.8926e-04
[c, p] = corr(stdMousePeakPredictionsMultiFiltTestGliresLiver, distanceFromMouseGlires, 'type', 'Spearman') % c = 0.1811, p = 0.1815
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiLiverStd=sortrows([distanceFromMouseGlires stdMousePeakPredictionsMultiFiltTestGliresLiver]);
sortedMultiLiverStd_x = sortedMultiLiverStd(:,1);
sortedMultiLiverStd_y = sortedMultiLiverStd(:,2);
yy1MultiStdLiver = smooth(sortedMultiLiverStd_x, sortedMultiLiverStd_y, 0.9, 'loess');
plot(sortedMultiLiverStd_x, sortedMultiLiverStd_y, 'k.', sortedMultiLiverStd_x, yy1MultiStdLiver, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for PV mouse-only model
mousePeakPredictionsInfoPV = importdata('print_speciesInfo.txt');
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,7)) % c = -0.8079, p = 5.2710e-14
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,7), 'type', 'Spearman') % c = -0.4585, p = 3.7946e-04
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMousePV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,7)]);
sortedMousePV_x = sortedMousePV(:,1);
sortedMousePV_y = sortedMousePV(:,2);
yy1MousePV = smooth(sortedMousePV_x, sortedMousePV_y, 0.9, 'loess');
plot(sortedMousePV_x, sortedMousePV_y, 'k.', sortedMousePV_x, yy1MousePV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.15738824, 0.15738824], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,9)) % c = 0.8743, p = 1.3722e-18
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,9), 'type', 'Spearman') % c = 0.8401, p = 5.7755e-16
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMouseStdPV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,9)]);
sortedMouseStdPV_x = sortedMouseStdPV(:,1);
sortedMouseStdPV_y = sortedMouseStdPV(:,2);
yy1MouseStdPV = smooth(sortedMouseStdPV_x, sortedMouseStdPV_y, 0.9, 'loess');
plot(sortedMouseStdPV_x, sortedMouseStdPV_y, 'k.', sortedMouseStdPV_x, yy1MouseStdPV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for PV mouse-only model with N
% filtering
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,8)) % c = -0.7887, p = 5.2750e-13
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,8), 'type', 'Spearman') % c = -0.4454, p = 5.8287e-04
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMousePV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,8)]);
sortedMousePV_x = sortedMousePV(:,1);
sortedMousePV_y = sortedMousePV(:,2);
yy1MousePV = smooth(sortedMousePV_x, sortedMousePV_y, 0.9, 'loess');
plot(sortedMousePV_x, sortedMousePV_y, 'k.', sortedMousePV_x, yy1MousePV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.15738824, 0.15738824], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,10)) % c = 0.8989, p = 5.3023e-21
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,10), 'type', 'Spearman') % c = 0.8726, p = 1.9466e-18
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMouseStdPV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,10)]);
sortedMouseStdPV_x = sortedMouseStdPV(:,1);
sortedMouseStdPV_y = sortedMouseStdPV(:,2);
yy1MouseStdPV = smooth(sortedMouseStdPV_x, sortedMouseStdPV_y, 0.9, 'loess');
plot(sortedMouseStdPV_x, sortedMouseStdPV_y, 'k.', sortedMouseStdPV_x, yy1MouseStdPV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for PV multi-species model
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,3)) % c = -0.8798, p = 4.3814e-19
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,3), 'type', 'Spearman') % c = -0.5124, p = 5.4092e-05
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMousePV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,3)]);
sortedMousePV_x = sortedMousePV(:,1);
sortedMousePV_y = sortedMousePV(:,2);
yy1MousePV = smooth(sortedMousePV_x, sortedMousePV_y, 0.9, 'loess');
plot(sortedMousePV_x, sortedMousePV_y, 'k.', sortedMousePV_x, yy1MousePV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.20393197, 0.20393197], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,6)) % c = 0.8848, p = 1.5010e-19
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,6), 'type', 'Spearman') % c = 0.6203, p = 3.4152e-07
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiStdPV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,6)]);
sortedMultiStdPV_x = sortedMultiStdPV(:,1);
sortedMultiStdPV_y = sortedMultiStdPV(:,2);
yy1MultiStdPV = smooth(sortedMultiStdPV_x, sortedMultiStdPV_y, 0.9, 'loess');
plot(sortedMultiStdPV_x, sortedMultiStdPV_y, 'k.', sortedMultiStdPV_x, yy1MultiStdPV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for PV multi-species model 
% with N filtering
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,4)) % c = -0.8727, p = 1.8822e-18
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,4), 'type', 'Spearman') % c = -0.5092, p = 6.1265e-05
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiPV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,4)]);
sortedMultiPV_x = sortedMultiPV(:,1);
sortedMultiPV_y = sortedMultiPV(:,2);
yy1MultiPV = smooth(sortedMultiPV_x, sortedMultiPV_y, 0.9, 'loess');
plot(sortedMultiPV_x, sortedMultiPV_y, 'k.', sortedMultiPV_x, yy1MultiPV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.20393197, 0.20393197], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,5)) % c = 0.9102, p = 2.5390e-22
[c, p] = corr(mousePeakPredictionsInfoPV.data(:,2), mousePeakPredictionsInfoPV.data(:,5), 'type', 'Spearman') % c = 0.6402, p = 1.0806e-07
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiStdPV=sortrows([mousePeakPredictionsInfoPV.data(:,2) mousePeakPredictionsInfoPV.data(:,5)]);
sortedMultiStdPV_x = sortedMultiStdPV(:,1);
sortedMultiStdPV_y = sortedMultiStdPV(:,2);
yy1MultiStdPV = smooth(sortedMultiStdPV_x, sortedMultiStdPV_y, 0.9, 'loess');
plot(sortedMultiStdPV_x, sortedMultiStdPV_y, 'k.', sortedMultiStdPV_x, yy1MultiStdPV, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for retina mouse-only model 
mousePeakPredictionsTestGliresRetina = importdata('mouse_only_scores_nofilter.txt', '\t', 1);
[c, p] = corr(mousePeakPredictionsTestGliresRetina.data(:,1), mousePeakPredictionsTestGliresRetina.data(:,2)) % c = -0.6975, p = 2.3188e-09
[c, p] = corr(mousePeakPredictionsTestGliresRetina.data(:,1), mousePeakPredictionsTestGliresRetina.data(:,2), 'type', 'Spearman') % c = -0.4545, p = 4.3349e-04
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedRetina=sortrows([mousePeakPredictionsTestGliresRetina.data(:,1) mousePeakPredictionsTestGliresRetina.data(:,2)]);
sortedRetina_x = sortedRetina(:,1);
sortedRetina_y = sortedRetina(:,2);
yy1Retina = smooth(sortedRetina_x, sortedRetina_y, 0.9, 'loess');
plot(sortedRetina_x, sortedRetina_y, 'k.', sortedRetina_x, yy1Retina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.14055635, 0.14055635], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
mousePeakStdsTestGliresRetina = importdata('mouse_only_stds_nofilter.txt', '\t', 1);
[c, p] = corr(mousePeakStdsTestGliresRetina.data(:,1), mousePeakStdsTestGliresRetina.data(:,2)) % c = 0.7282, p = 2.0224e-10
[c, p] = corr(mousePeakStdsTestGliresRetina.data(:,1), mousePeakStdsTestGliresRetina.data(:,2), 'type', 'Spearman') % c = 0.4331, p = 8.5604e-04
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedStdRetina=sortrows([mousePeakStdsTestGliresRetina.data(:,1) mousePeakStdsTestGliresRetina.data(:,2)]);
sortedStdRetina_x = sortedStdRetina(:,1);
sortedStdRetina_y = sortedStdRetina(:,2);
yy1StdRetina = smooth(sortedStdRetina_x, sortedStdRetina_y, 0.9, 'loess');
plot(sortedStdRetina_x, sortedStdRetina_y, 'k.', sortedStdRetina_x, yy1StdRetina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for retina mouse-only model 
% with N filtering
mousePeakPredictionsTestGliresRetina = importdata('mouse_only_scores_Nfiltered.txt', '\t', 1);
[c, p] = corr(mousePeakPredictionsTestGliresRetina.data(:,1), mousePeakPredictionsTestGliresRetina.data(:,2)) % c = -0.7507, p = 2.7050e-11
[c, p] = corr(mousePeakPredictionsTestGliresRetina.data(:,1), mousePeakPredictionsTestGliresRetina.data(:,2), 'type', 'Spearman') % c = -0.4668, p = 2.8692e-04
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedRetina=sortrows([mousePeakPredictionsTestGliresRetina.data(:,1) mousePeakPredictionsTestGliresRetina.data(:,2)]);
sortedRetina_x = sortedRetina(:,1);
sortedRetina_y = sortedRetina(:,2);
yy1Retina = smooth(sortedRetina_x, sortedRetina_y, 0.9, 'loess');
plot(sortedRetina_x, sortedRetina_y, 'k.', sortedRetina_x, yy1Retina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.14055635, 0.14055635], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
mousePeakStdsTestGliresRetina = importdata('mouse_only_stds_Nfiltered.txt', '\t', 1);
[c, p] = corr(mousePeakStdsTestGliresRetina.data(:,1), mousePeakStdsTestGliresRetina.data(:,2)) % c = 0.7328, p = 1.3600e-10
[c, p] = corr(mousePeakStdsTestGliresRetina.data(:,1), mousePeakStdsTestGliresRetina.data(:,2), 'type', 'Spearman') % c = 0.4412, p = 6.6605e-04
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedStdRetina=sortrows([mousePeakStdsTestGliresRetina.data(:,1) mousePeakStdsTestGliresRetina.data(:,2)]);
sortedStdRetina_x = sortedStdRetina(:,1);
sortedStdRetina_y = sortedStdRetina(:,2);
yy1StdRetina = smooth(sortedStdRetina_x, sortedStdRetina_y, 0.9, 'loess');
plot(sortedStdRetina_x, sortedStdRetina_y, 'k.', sortedStdRetina_x, yy1StdRetina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for retina multi-species model 
mousePeakPredictionsMultiTestGliresRetina = importdata('multispecies_scores_nofilter.txt', '\t', 1);
[c, p] = corr(mousePeakPredictionsMultiTestGliresRetina.data(:,1), mousePeakPredictionsMultiTestGliresRetina.data(:,2)) % c = -0.6381, p = 1.2219e-07
[c, p] = corr(mousePeakPredictionsMultiTestGliresRetina.data(:,1), mousePeakPredictionsMultiTestGliresRetina.data(:,2), 'type', 'Spearman') % c = -0.4064, p = 0.0019
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiRetina=sortrows([mousePeakPredictionsMultiTestGliresRetina.data(:,1) mousePeakPredictionsMultiTestGliresRetina.data(:,2)]);
sortedMultiRetina_x = sortedMultiRetina(:,1);
sortedMultiRetina_y = sortedMultiRetina(:,2);
yy1MultiRetina = smooth(sortedMultiRetina_x, sortedMultiRetina_y, 0.9, 'loess');
plot(sortedMultiRetina_x, sortedMultiRetina_y, 'k.', sortedMultiRetina_x, yy1MultiRetina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.12042443, 0.12042443], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
mousePeakStdsMultiTestGliresRetina = importdata('multispecies_stds_nofilter.txt', '\t', 1);
[c, p] = corr(mousePeakStdsMultiTestGliresRetina.data(:,1), mousePeakStdsMultiTestGliresRetina.data(:,2)) % c = 0.7434, p = 5.2776e-11
[c, p] = corr(mousePeakStdsMultiTestGliresRetina.data(:,1), mousePeakStdsMultiTestGliresRetina.data(:,2), 'type', 'Spearman') % c = 0.4555, p = 4.1876e-04
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiStdRetina=sortrows([mousePeakStdsMultiTestGliresRetina.data(:,1) mousePeakStdsMultiTestGliresRetina.data(:,2)]);
sortedMultiStdRetina_x = sortedMultiStdRetina(:,1);
sortedMultiStdRetina_y = sortedMultiStdRetina(:,2);
yy1MultiStdRetina = smooth(sortedMultiStdRetina_x, sortedMultiStdRetina_y, 0.9, 'loess');
plot(sortedMultiStdRetina_x, sortedMultiStdRetina_y, 'k.', sortedMultiStdRetina_x, yy1MultiStdRetina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off

% Plot evolutionary distance vs. predictions for retina multi-species model 
% with N filtering
mousePeakPredictionsMultiTestGliresRetina = importdata('multispecies_scores_Nfiltered.txt', '\t', 1);
[c, p] = corr(mousePeakPredictionsMultiTestGliresRetina.data(:,1), mousePeakPredictionsMultiTestGliresRetina.data(:,2)) % c = -0.6697, p = 1.6628e-08
[c, p] = corr(mousePeakPredictionsMultiTestGliresRetina.data(:,1), mousePeakPredictionsMultiTestGliresRetina.data(:,2), 'type', 'Spearman') % c = -0.4129, p = 0.0016
xlim([0 82])
hold on
ylim([0 1])
hold on
sortedMultiRetina=sortrows([mousePeakPredictionsMultiTestGliresRetina.data(:,1) mousePeakPredictionsMultiTestGliresRetina.data(:,2)]);
sortedMultiRetina_x = sortedMultiRetina(:,1);
sortedMultiRetina_y = sortedMultiRetina(:,2);
yy1MultiRetina = smooth(sortedMultiRetina_x, sortedMultiRetina_y, 0.9, 'loess');
plot(sortedMultiRetina_x, sortedMultiRetina_y, 'k.', sortedMultiRetina_x, yy1MultiRetina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold on
line([0, 82], [0.12042443, 0.12042443], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
hold off
mousePeakStdsMultiTestGliresRetina = importdata('multispecies_stds_Nfiltered.txt', '\t', 1);
[c, p] = corr(mousePeakStdsMultiTestGliresRetina.data(:,1), mousePeakStdsMultiTestGliresRetina.data(:,2)) % c = 0.7411, p = 6.5508e-11
[c, p] = corr(mousePeakStdsMultiTestGliresRetina.data(:,1), mousePeakStdsMultiTestGliresRetina.data(:,2), 'type', 'Spearman') % c = 0.4525, p = 4.6330e-04
xlim([0 82])
hold on
ylim([0.2 0.4])
hold on
sortedMultiStdRetina=sortrows([mousePeakStdsMultiTestGliresRetina.data(:,1) mousePeakStdsMultiTestGliresRetina.data(:,2)]);
sortedMultiStdRetina_x = sortedMultiStdRetina(:,1);
sortedMultiStdRetina_y = sortedMultiStdRetina(:,2);
yy1MultiStdRetina = smooth(sortedMultiStdRetina_x, sortedMultiStdRetina_y, 0.9, 'loess');
plot(sortedMultiStdRetina_x, sortedMultiStdRetina_y, 'k.', sortedMultiStdRetina_x, yy1MultiStdRetina, 'k-', 'MarkerSize', 11, 'LineWidth', 1)
hold off