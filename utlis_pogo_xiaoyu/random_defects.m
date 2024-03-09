classdef random_defects
    properties
        numCircles
        circles % Each row: [x, y, radius, depth, thickness]
        numPoints = 100 % Points to define circumference
    end
    
    methods
        % Constructor
        function obj = random_defects(numCircles)
            obj.numCircles = numCircles;
            obj.circles = zeros(numCircles, 5); % Initialize circles array
        end

        % Method to generate circles
        function obj = generateCircles(obj)
            for i = 1:obj.numCircles
                overlap = true;
                while overlap
                    x = 0.5 + 5 * rand(); % X-coordinate
                    y = 0.5 + 5 * rand(); % Y-coordinate
                    radius = 0.5 + 2 * rand(); % Radius
                    % depth = randi([1, 7]); % Depth
					depth = 1 + 6 * rand(); % Depth
                    thickness = 10 + 10 * rand();  % Thickness

                    newCircle = [x, y, radius, depth, thickness];
                    overlap = false;

                    % Check for overlap
                    for j = 1:i-1
                        if obj.isOverlapping(newCircle, obj.circles(j, :))
                            overlap = true;
                            break;
                        end
                    end
                end
                obj.circles(i, :) = newCircle; % Add non-overlapping circle
            end
        end

        % Method to plot circles
        function plotCircles(obj)
            figure;
            hold on;
            axis equal;
            xlim([0, 6]);
            ylim([0, 6]);
            theta = linspace(0, 2*pi, obj.numPoints);

            for i = 1:obj.numCircles
                xCenter = obj.circles(i, 1);
                yCenter = obj.circles(i, 2);
                radius = obj.circles(i, 3);

                % Circle points
                x = radius * cos(theta) + xCenter;
                y = radius * sin(theta) + yCenter;

                % Draw the circle
                plot(x, y, 'b');
            end

            hold off;
        end
    end

    methods (Access = private)
        % Helper method to check overlap
        function overlap = isOverlapping(~, circle1, circle2)
            distanceCenters = sqrt((circle1(1) - circle2(1))^2 + (circle1(2) - circle2(2))^2);
            sumRadii = circle1(3) + circle2(3);
            overlap = distanceCenters < sumRadii;
        end
    end
end
