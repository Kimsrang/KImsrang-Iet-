function [clusters, centroids] = k_means_cluster(k, points, maxIter)
    % Custom K-means clustering algorithm
    if nargin < 3
        maxIter = 100;
    end

    n = size(points, 1);            % number of points
    rng('default');                 % for reproducibility
    rand_idx = randperm(n, k);      % random initial centroids
    centroids = points(rand_idx, :);

    converged = false;
    iter = 0;

    while ~converged && iter < maxIter
        iter = iter + 1;
        clusters = cell(1, k);      % initialize cluster storage

        % Assign points to the nearest centroid
        for i = 1:n
            dists = vecnorm(points(i,:) - centroids, 2, 2);
            [~, cid] = min(dists);
            clusters{cid}(end+1) = i;  % store index of point
        end

        % Recalculate centroids
        new_centroids = zeros(k, 2);
        for j = 1:k
            if isempty(clusters{j})
                new_centroids(j, :) = points(randi(n), :);  % re-init empty cluster
            else
                new_centroids(j, :) = mean(points(clusters{j}, :), 1);
            end
        end

        % Check convergence
        if isequal(new_centroids, centroids)
            converged = true;
        end

        centroids = new_centroids;
    end
end
