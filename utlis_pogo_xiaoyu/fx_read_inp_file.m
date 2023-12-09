function data = fx_read_inp_file(file_name)
% Open the file
fid = fopen(file_name, 'rt');
if fid == -1
    error('Could not open the file.');
end

%     % Define global variable
%     global currentLine;
%     global totalLines;

%     currentLine = 0;
%     % Get the total lines for progress calculation
%     frewind(fid);
%     totalLines = numel(cell2mat(textscan(fid,'%1c%*[^\n]'))); % count '\n' occurrences
%     frewind(fid); % rewind the file after the line count

% Initialize the waitbar
%     h = waitbar(0, 'Please wait... Reading file.');

% Initialize data
data = struct();

while ~feof(fid)
    line = fgetl(fid);
    %         currentLine = currentLine + 1;
    %         % Update the waitbar
    %         waitbar(currentLine / totalLines, h, sprintf('Line %d out of %d', currentLine, totalLines));
    % ************** full reading, comment out for fast read *************
    if contains(line, '*Node')
        disp("Node");
        % data.nodes = read_nodes(fid);
        data.nodes = read_nodes_maxmin(fid);
    elseif contains(line, '*Element')
        disp("Element");
        % data.elements = read_elements(fid);
    elseif contains(line, '*ElSet')
        disp("Elset");
        elset_name = regexp(line, 'ElSet=(\w+)', 'tokens');
        data.elsets.(elset_name{1}{1}) = read_elset(fid);
    elseif contains(line, '*Material')
        break; % stop reading after this ...
        material_name = regexp(line, 'Name=(\w+)', 'tokens');
        data.materials.(material_name{1}{1}) = read_material(fid);
    end
    % ******************** only read elements
    % if contains(line, '*Node')
    %     disp("Node");
    %     continue;
    % elseif contains(line, '*Element')
    %     disp("Element");
    %     continue;
    % elseif contains(line, '*ElSet')
    %     elset_name = regexp(line, 'ElSet=(\w+)', 'tokens');
    %     data.elsets.(elset_name{1}{1}) = read_elset(fid);
    % elseif contains(line, '*Material')
    %     break; % stop immediately
    % end
end
%     % Close the waitbar
%     close(h);

% Close the file
fclose(fid);

    function nodes = read_nodes(fid)
        nodes = NaN(1e5, 3);  % preallocate with a large size
        n = 0;  % number of nodes so far

        while ~feof(fid)
            prevPos = ftell(fid);
            line = fgetl(fid);
            %             currentLine = currentLine + 1;
            %             % Update the waitbar
            %             waitbar(currentLine / totalLines, h, sprintf('Line %d out of %d', currentLine, totalLines));

            if startsWith(line, '*')
                fseek(fid, prevPos, 'bof');  % move the file pointer back to the start of the line
                break;
            end

            node = sscanf(line, '%*d, %f, %f, %f');
            n = n + 1;
            if n > size(nodes, 1)  % need to grow the array
                nodes = [nodes; NaN(1e5, 3)];  % grow by a large chunk
            end
            nodes(n, :) = node';
        end

        nodes = nodes(1:n, :);  % remove extra preallocated space
    end

    % *********************
    function nodes = read_nodes_maxmin(fid)
        nodes_max = ones(3, 1)*-1e10; 
        nodes_min = ones(3, 1)*1e10; 
        while ~feof(fid)
            prevPos = ftell(fid);
            line = fgetl(fid);
            %             currentLine = currentLine + 1;
            %             % Update the waitbar
            %             waitbar(currentLine / totalLines, h, sprintf('Line %d out of %d', currentLine, totalLines));
            if startsWith(line, '*')
                fseek(fid, prevPos, 'bof');  % move the file pointer back to the start of the line
                break;
            end

            node = sscanf(line, '%*d, %f, %f, %f');
            nodes_max = max(nodes_max, node);
            nodes_min = min(nodes_min, node);
        end

        nodes = [nodes_max nodes_min];  % remove extra preallocated space
    end
    % *********************

    function elements = read_elements(fid)
        elements = NaN(1e5, 8);
        n = 0;

        while ~feof(fid)
            prevPos = ftell(fid);
            line = fgetl(fid);
            %             currentLine = currentLine + 1;
            %             % Update the waitbar
            %             waitbar(currentLine / totalLines, h, sprintf('Line %d out of %d', currentLine, totalLines));

            if startsWith(line, '*')
                fseek(fid, prevPos, 'bof');  % move the file pointer back to the start of the line
                break;
            end

            element = sscanf(line, '%*d, %d, %d, %d, %d, %d, %d, %d, %d'); % default 8 nodes
            n = n + 1;
            if n > size(elements, 1)
                elements = [elements; NaN(1e5, 8)];
            end
            elements(n, :) = element';
        end

        elements = elements(1:n, :);
    end

    function elset = read_elset(fid)
        elset = int32(zeros(1, 1e6));
        % elset = nan(1, 1e6);
        n = 0;

        while ~feof(fid)
            prevPos = ftell(fid);
            line = fgetl(fid);
            %             currentLine = currentLine + 1;
            %             % Update the waitbar
            %             waitbar(currentLine / totalLines, h, sprintf('Line %d out of %d', currentLine, totalLines));
            if startsWith(line, '*')
                fseek(fid, prevPos, 'bof');  % move the file pointer back to the start of the line
                break;
            end
            elset_data = sscanf(line, '%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d');
            dn = length(elset_data);
            n = n + dn;
            if n > size(elset, 2)
                elset = [elset int32(zeros(1, 1e6))];
                disp('more than 1e6 reached')
            end
            elset(1, n-dn+1:n) = elset_data;
        end

        elset = elset(1, 1:n);
    end

    function material = read_material(fid)
        material = NaN(1, 1e4);
        n = 0;

        while ~feof(fid)
            prevPos = ftell(fid);
            line = fgetl(fid);
            %             currentLine = currentLine + 1;
            %             % Update the waitbar
            %             waitbar(currentLine / totalLines, h, sprintf('Line %d out of %d', currentLine, totalLines));

            if startsWith(line, '*')
                fseek(fid, prevPos, 'bof');  % move the file pointer back to the start of the line
                break;
            end

            material_data = sscanf(line, '%f, %f, %f, %f, %f, %f, %f, %f, %f');
            n = n + 1;
            if n > size(material, 2)
                material = [material NaN(1, 1e4)];
            end
            material(1, n) = material_data;
        end

        material = material(1, 1:n);
    end
end
