function fx_display_model(model, flag_plotelements)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(model.matTypeRefs)
    disp('matTypeRefs is empty');
else
    figure,
    plot(model.matTypeRefs(1:10:end));
end

if flag_plotelements
    nodes    = model.nodePos.';
    elements = model.elNodes.';

    figure;
    dims = size(nodes, 2);
    unique_mat = unique(model.matTypeRefs);
    % check if 2D
    if dims==2
        for i = 2:length(unique_mat)
            elsets = find(model.matTypeRefs==unique_mat(i)).';
            centroids = zeros(size(elsets, 2), 2);  % Preallocate centroid array
            for i = 1:1:size(elsets, 2) % skip by ? to save time
                element_nodes   = elements(elsets(:, i), 1:end); % the first colume includeds
                nodes_pos       = nodes(element_nodes, 1:end);
                centroids(i, :) = mean(nodes_pos, 1);
            end
            % skip to save memory
            X = centroids(1:1:end, 1);
            Y = centroids(1:1:end, 2);
            % scatter3(X, Y, Z, 'filled', 'MarkerFaceAlpha', trans);
            plot(X, Y, '.', 'MarkerSize', 1);
            hold on;
        end
    % check 3D or higher
    else
        for i = 2:length(unique_mat)
            elsets = find(model.matTypeRefs==unique_mat(i)).';
            centroids = zeros(size(elsets, 2), 3);  % Preallocate centroid array
            for i = 1:10:size(elsets, 2) % skip by 10 to save time
                element_nodes   = elements(elsets(:, i), 1:end); % the first colume includeds
                nodes_pos       = nodes(element_nodes, 1:end);
                centroids(i, :) = mean(nodes_pos, 1);
            end
            % skip to save memory
            X = centroids(1:1:end, 1);
            Y = centroids(1:1:end, 2);
            Z = centroids(1:1:end, 3);

            % scatter3(X, Y, Z, 'filled', 'MarkerFaceAlpha', 0.5);
            plot3(X, Y, Z, '.', 'MarkerSize', 1);
            hold on;
        end
        % plot generator
        if isfield(model, 'shots')
            element_nodes = model.shots{1, 1}.sigs{1, 1}.nodeSpec;
            centroids     = nodes(element_nodes, 1:end);
            X = centroids(:, 1);
            Y = centroids(:, 2);
            Z = centroids(:, 3);
            plot3(X, Y, Z, 'o', 'MarkerSize', 1);
            hold on;
        end
        % plot receiver
        if isfield(model, 'measSets')
            element_nodes = model.measSets{1, 1}.measNodes;
            centroids     = nodes(element_nodes, 1:end);
            X = centroids(:, 1);
            Y = centroids(:, 2);
            Z = centroids(:, 3);
            plot3(X, Y, Z, 'o', 'MarkerSize', 1);
            hold on;
        end
        % plot fixed fixNodes
        if isfield(model, 'fixNodes')
            element_nodes = model.fixNodes;
            centroids     = nodes(element_nodes, 1:end);
            X = centroids(:, 1);
            Y = centroids(:, 2);
            Z = centroids(:, 3);
            plot3(X, Y, Z, 'o', 'MarkerSize', 1);
            hold on;
        end
    end

axis equal;  % make the axes scale equal
else
    disp('do not display elements');
end

end

