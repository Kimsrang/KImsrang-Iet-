clc; clear all;

%% --- Load Data ---
filename = 'Load Data1.xlsx';
data = readtable(filename);

% Limit to 129 rows max
data = data(1:min(129, height(data)), :);

% Extract relevant data
x = data.X;
y = data.Y;
points = [x, y];
loadValues = data.P_kW__Pf_0_95;
phaseType = data.Phase;
numLoads = length(loadValues);

%% Phase Arrangement for Each Load
[phases, phaseCurrents, deviations, UC] = load_phase_balancing(filename);


%% --- K-means Clustering ---
k = 15;
[clusters, centroids] = k_means_cluster(k, points);

% Remove NaNs from input data
P_1ph = data.P_kW__Pf_0_95;
valid = ~(isnan(x) | isnan(y) | isnan(P_1ph));
x = x(valid);
y = y(valid);
P_1ph = P_1ph(valid);
points = [x(1:129), y(1:129)];

% %% --- TRANSFORMER LOCATION ---
x =100;
y=400;
best_coord = [x,y];

%% Define line segments for pole generation
segments = {
    [100, 500; 400, 500];
    [100.75, 425.75; 100, 500];
    [48.87, 662; 2.1, 28.79];
    [126.8, 400.8; 230, 400];
    [2.1, 28.79; 480, 28.79];
    [77.25, 401.66; 29.75, 403]
};

max_spacing = 30;
tolerance = 1e-3;
min_dist = 1e-2;

all_poles = [];

for i = 1:length(segments)
    start_point = segments{i}(1, :);
    end_point   = segments{i}(2, :);

    seg_length = norm(end_point - start_point);
    num_segments = ceil(seg_length / max_spacing);
    t = linspace(0, 1, num_segments + 1);
    poles = (1 - t') * start_point + t' * end_point;
    all_poles = [all_poles; poles];
end

% Remove duplicates
all_poles_rounded = round(all_poles / tolerance) * tolerance;
[unique_poles, ~, ~] = unique(all_poles_rounded, 'rows');

% Remove poles too close to transformer or loads
dist_to_transformer = sqrt((unique_poles(:,1) - best_coord(1)).^2 + (unique_poles(:,2) - best_coord(2)).^2);
min_dist_to_load = zeros(size(unique_poles,1),1);
for i = 1:size(unique_poles,1)
    dx = x - unique_poles(i,1);
    dy = y - unique_poles(i,2);
    d = sqrt(dx.^2 + dy.^2);
    min_dist_to_load(i) = min(d);
end
valid_idx = (dist_to_transformer > min_dist) & (min_dist_to_load > min_dist);
filtered_poles = unique_poles(valid_idx, :);


%% --- Plot: Clusters, Phases, Boundaries, and Segments ---
figure; hold on;
cluster_colors = lines(k);
markerSize = 60;

% Phase color map
phaseColor = {
    [1 0 0];  % A - red
    [0 1 0];  % B - green
    [0 0 1];  % C - blue
    [0 0 0];  % 3-phase - black
};

% Load connection to nearest pole (Freely Strategy)
nearest_idx = zeros(numLoads, 1);
for i = 1:numLoads
    distances = sqrt(sum((unique_poles - points(i, :)).^2, 2));
    [~, nearest_idx(i)] = min(distances);
end

for i = 1:numLoads
    if loadValues(i) == 0
        continue;  % Skip loads with zero power
    end
    pt = points(i, :);                         % Load coordinate
    np = unique_poles(nearest_idx(i), :);     % Nearest pole

    % Determine phase color
    ph = phases(i);
    if ph == 0
        c = phaseColor{4};
    else
        c = phaseColor{ph};
    end

    % Plot connection
    plot([pt(1), np(1)], [pt(2), np(2)], '-', 'Color', c, 'LineWidth', 1.5);
end

%----------------------------------------------------------------------
%%% put lable LV poles
filename1 = 'Feeders.xlsx';
sheetname2 = 'Shortest-Path Strategy';

% Read data
data1 = readtable(filename1, 'Sheet', sheetname2);
data1 = data1(1:min(47, height(data1)), :);  % fixed variable name

% Extract columns
x1 = data1.X1;
y1 = data1.Y1;
ith_pole = data1.ithLVPole;
ith_feeder = data1.ithFeeder;

% Unique feeders and color map
unique_feeders = unique(ith_feeder);
colors = lines(length(unique_feeders));

for n = 1:length(unique_feeders)
    feeder_id = unique_feeders(n);
    idx = ith_feeder == feeder_id;

    x_group = x1(idx);
    y_group = y1(idx);
    poles = ith_pole(idx);

    % Place labels only (no points)
    labels = string(poles);
    label_color = [1 0 0];
    text(x_group, y_group, labels, ...
         'FontSize', 12, ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top', ...
         'Color', label_color);
end

% Draw feeder boundary
feeders = {
    [10, 45, 70, 30, 10], [415, 415, 670, 670, 418], -30, 600, 'Feeder-04';
    [85, 10, -25, 200, 200, 30, 60, 85, 85], [410, 410, 5, 5, 45, 45, 380, 380, 410], 100, 100, 'Feeder-03';
    [110, 90, 90, 420, 420, 110, 110], [410, 410, 520, 520, 480, 480, 410], 300, 470, 'Feeder-02';
    [120, 240, 240, 120, 120], [380, 380, 410, 410, 380], 200, 360, 'Feeder-01'
};

for i = 1:size(feeders,1)
    x = feeders{i,1}; y = feeders{i,2};
    fill(x, y, [0 0 1], 'FaceAlpha', 0.08, 'EdgeColor', 'none')
    plot(x, y, 'k--', 'LineWidth', 1.5)
    text(feeders{i,3}, feeders{i,4}, feeders{i,5}, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', 'k')
end

% ---------------------------------------------------------------------
% --- Plot each cluster ---
for i = 1:k
    idx = clusters{i};             % indices in this cluster
    cluster_pts = points(idx, :);  % points in this cluster
    cluster_phases = phases(idx);  % phases for those points

    % Plot each point with phase color
    for j = 1:length(idx)
        ph = cluster_phases(j);
        if ph == 0
            c = phaseColor{4};
        else
            c = phaseColor{ph};
        end
        scatter(cluster_pts(j,1), cluster_pts(j,2), markerSize, 'filled', ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    end

end


% --- Optional: Plot line segments (e.g., network links) ---
% segments should be a cell array of [x1 y1; x2 y2] format
if exist('segments', 'var')
    for i = 1:length(segments)
        seg = segments{i};
        plot(seg(:,1), seg(:,2), '-', 'Color', [0 0 1], 'LineWidth', 1); % Blue lines
    end
end

%% Plot poles
scatter(unique_poles(:,1), unique_poles(:,2), 80, 'sb', 'filled');

%% --- Connect Key Poles to Transformer Location ---
% Coordinates of the specific poles
keyPoles = [
    100.75, 425.75;
     77.25, 401.66;
    126.80, 400.80
];

% Plot lines from each key pole to the transformer
for i = 1:size(keyPoles, 1)
    plot([keyPoles(i,1), best_coord(1)], ...
         [keyPoles(i,2), best_coord(2)], ...
         'b-', 'LineWidth', 1);  % magenta line to transformer
end

% Optional: mark these special poles distinctly
hKeyPoles = scatter(keyPoles(:,1), keyPoles(:,2), ...
                    120, 'md', 'filled', 'DisplayName', 'Key Poles');

% Plot transformer location
hTx = plot(100,400, 'kp', 'MarkerSize', 20, 'MarkerFaceColor', 'g', 'DisplayName', 'Optimal Transformer');

% Add transformer image
[img, ~, alpha] = imread('Transformer.png');
img = flipud(img); alpha = flipud(alpha);
scale = 80; offset_x = 170; offset_y = 250;
image('CData', img, 'XData', [offset_x, offset_x + scale], ...
      'YData', [offset_y, offset_y + scale], 'AlphaData', alpha);

% Arrow from image to transformer
image_center_x = offset_x + scale/2;
image_center_y = offset_y + scale/2;
quiver(image_center_x, image_center_y, best_coord(1)-image_center_x, best_coord(2)-image_center_y, ...
       0, 'b', 'LineWidth', 1, 'MaxHeadSize', 1);

%--------------------------------------------------------------------------
% Axis labels and title
xlabel('Longitude coordinate (X)');
ylabel('Latitude coordinate (Y)');
title('LV Distribution Network');

% Dummy scatter plots to create legend
hA = scatter(NaN, NaN, markerSize, 'filled', 'MarkerFaceColor', phaseColor{1}, 'MarkerEdgeColor', 'k');
hB = scatter(NaN, NaN, markerSize, 'filled', 'MarkerFaceColor', phaseColor{2}, 'MarkerEdgeColor', 'k');
hC = scatter(NaN, NaN, markerSize, 'filled', 'MarkerFaceColor', phaseColor{3}, 'MarkerEdgeColor', 'k');
ho = scatter(NaN, NaN, markerSize, 'filled', 'MarkerFaceColor', phaseColor{4}, 'MarkerEdgeColor', 'k');
hPole = scatter(unique_poles(:,1), unique_poles(:,2), 80, 'sb', 'filled');

legend([hA, hB, hC, ho, hPole, hKeyPoles, hTx], {'Phase A', 'Phase B', 'Phase C', '3phase-Load','Poles','Pole-Transformer', 'Transformer'}, 'Location', 'best');

grid on;
axis equal;
hold off;




