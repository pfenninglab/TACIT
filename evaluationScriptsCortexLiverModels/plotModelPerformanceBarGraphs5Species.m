% Plot test set AUCs, AUNPV-Specs., and AUPRCs for motor cortex mouse-only 
% model, lineage-specific OCR evaluations
c = colormap('jet');
X = categorical({'AUC','AUPRC', 'AUNPV-Spec.'});
ylim([0 1])
hold on
b = bar(X, [0.899807, 0.832461, 0.900778, 0.806337, 0.83813, 0.727621; 0.914656, 0.878085, 0.950812, 0.908861, 0.704652, 0.653692; 0.822464, 0.766648, 0.807049, 0.635652, 0.849879, 0.910827]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .1250, 1];
b(3).FaceColor = [0, .6875, 1];
b(4).FaceColor = [.2500, 1, 0.7500];
b(5).FaceColor = [.8125, 1, 0.1875];
b(6).FaceColor = [1, .6250, 0];
hold off

% Plot test set AUCs, AUNPV-Specs., and AUPRCs for PV mouse-only model, 
% lineage-specific OCR evaluations
ylim([0 1])
hold on
b = bar(X, [0.8633256568, 0.7645045997, 0.7232776355, 0.8869847462; 0.6337729431, 0.3799061124, 0.9901980793, 0.289803291; 0.9622170498, 0.9486087137, 0.07258715073, 0.9968230893]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .1250, 1];
b(3).FaceColor = [0, .6875, 1];
b(4).FaceColor = [.2500, 1, 0.7500];
hold off

% Plot test set AUCs and AUPRCs for retina mouse-only model,
% lineage-specific OCR evaluations
ylim([0 1])
hold on
b = bar(X,[0.939809201, 0.9337827258, 0.8469909962, 0.6297871273; 0.8361684185, 0.8304316956, 0.8017783946, 0.6164805372; 0.9645830968, 0.9621628777, 0.8354838226, 0.6976156669]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .1250, 1];
b(3).FaceColor = [0, .6875, 1];
b(4).FaceColor = [.2500, 1, 0.7500];
hold off

% Plot test set AUCs, AUNPV-Specs., and AUPRCs (listed last) for motor 
% cortex mouse-only model, tissue-specific OCR evaluations
xlim([0.5 3.5])
hold on
ylim([0.3 1])
hold on
scatter([1, 2, 3], [0.915813, 0.959043, 0.834725], 200, [0.50196078431, 0.50196078431, 0], '*', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.925546, 0.959277, 0.862012], 200, [0.50196078431, 0.50196078431, 0], '.', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.884518, 0.888927, 0.878646], 200, [0.50196078431, 0.50196078431, 0], 'x', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.775727, 0.781207, 0.746957], 200, [0.5000, 0, 0], '*', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.763, 0.650992, 0.83417], 200, [0.5000, 0, 0], '.', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.674043, 0.381469, 0.854495], 200, [0.5000, 0, 0], 'x', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.788062, 0.753016, 0.801051], 200, [0.5000, 0, 0], '+', 'jitter', 'on', 'jitterAmount', 0.24)
hold off

% Plot test set AUCs and AUPRCs for PV mouse-only model, tissue-specific 
% OCR evaluations
xlim([0.5 3.5])
hold on
ylim([0.3 1])
hold on
scatter([1, 2, 3], [0.8717390387, 0.90764975147924, 0.8312851897], 200, [0.50196078431, 0.50196078431, 0], '*', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.8080962364, 0.800407743956581, 0.8035023687], 200, [0.5000, 0, 0], '*', 'jitter', 'on', 'jitterAmount', 0.24)
hold on
scatter([1, 2, 3], [0.7838345865, 0.496364097, 0.924803131700812], 200, [0.5000, 0, 0], '.', 'jitter', 'on', 'jitterAmount', 0.24)
hold off

% Plot test set AUCs and AUPRCs for retina mouse-only model, 
% tissue-specific OCR evaluations
xlim([0.5 3.5])
hold on
ylim([0.3 1])
hold on
scatter([1, 2, 3], [0.8608187483, 0.898736068, 0.7945431055], 200, [0.50196078431, 0.50196078431, 0], '*', 'jitter', 'on', 'jitterAmount', 0.24)
hold off

% Plot test set AUCs, AUNPV-Specs., and AUPRCs for motor cortex 
% multi-species model lineage-specific and tissue-specific OCR evaluations
c = colormap('jet');
X = categorical({'AUC','AUPRC', 'AUNPV-Spec.'});
ylim([0 1])
hold on
b = bar(X, [0.908534, 0.935442, 0.872285, 0.922831, 0.848658, 0.930891, 0.844142; 0.901379, 0.967784, 0.9339, 0.848603, 0.84214, 0.898375, 0.872848; 0.915066, 0.878459, 0.760258, 0.959772, 0.850448, 0.956803, 0.804923]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .6875, 1];
b(3).FaceColor = [.2500, 1, 0.7500];
b(4).FaceColor = [.8125, 1, 0.1875];
b(5).FaceColor = [1, .6250, 0];
b(6).FaceColor = [0.5000, 0.5000, 0];
b(7).FaceColor = [0.5000, 0, 0];
hold off

% Plot test set AUCs, AUNPV-Specs., and AUPRCs for liver multi-species 
% model, lineage-specific and tissue-specific OCR evaluations
c = colormap('jet');
X = categorical({'AUC','AUPRC', 'AUNPV-Spec.'});
ylim([0 1])
hold on
b = bar(X, [0.884308, 0.85779, 0.774324, 0.848686, 0.818934, 0.930461; 0.834676, 0.935833, 0.94332, 0.825578, 0.813744, 0.80999; 0.915824, 0.741417, 0.462808, 0.838509, 0.812046, 0.978801]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .6875, 1];
b(3).FaceColor = [.2500, 1, 0.7500];
b(4).FaceColor = [.8125, 1, 0.1875];
b(5).FaceColor = [1, .6250, 0];
b(6).FaceColor = [1, 0.0625, 0];
hold off

% Plot test set AUCs, AUNPV-Specs., and AUPRCs for PV multi-species model 
% lineage-specific and tissue-specific OCR evaluations
c = colormap('jet');
X = categorical({'AUC','AUPRC', 'AUNPV-Spec.'});
ylim([0 1])
hold on
b = bar(X, [0.8794664103, 0.7734464757, 0.91343153, 0.9125464345, 0.8455555563; 0.5541698054, 0.9923158698, 0.34601527, 0.8951816634, 0.846464314; 0.9795461036, 0.09697814239, 0.997637366, 0.9328201021, 0.8365010645]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .6875, 1];
b(3).FaceColor = [.2500, 1, 0.7500];
b(4).FaceColor = [0.5000, 0.5000, 0];
b(5).FaceColor = [0.5000, 0, 0];
hold off

% Plot test set AUCs and AUPRCs for retina multi-species model 
% lineage-specific and tissue-specific OCR evaluations
c = colormap('jet');
X = categorical({'AUC','AUPRC', 'AUNPV-Spec.'});
ylim([0 1])
hold on
b = bar(X, [0.9522765658, 0.8025796036, 0.7487094003, 0.8244066464; 0.8619156204, 0.745493866244806, 0.7148984873, 0.916803768494305; 0.963165064974524, 0.664744090853522, 0.766816122551917, 0.664744090853522]);
b(1).FaceColor = [0, 0, .5625];
b(2).FaceColor = [0, .6875, 1];
b(3).FaceColor = [.2500, 1, 0.7500];
b(4).FaceColor = [0.5000, 0.5000, 0];
hold off