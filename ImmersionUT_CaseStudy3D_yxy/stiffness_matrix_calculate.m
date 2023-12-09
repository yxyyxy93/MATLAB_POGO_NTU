% using struct and class
close all;
fclose all;
clear;
clc;

% % define and calculate the properties
elastic_constants_60per = [
    135.5e9 9.65e9 9.65e9 ...
    0.164 0.164 0.480 ...
    5.31e9 5.31e9 3.23e9
    ];

% from chamis model in Ruben's paper
elastic_constants = [
    223.366e9 18.107e9 18.107e9 ...
    0.006 0.006 0.599 ...
    51.951e9 51.951e9 5.661e9
    ];

Sample                = struct(...
    'c22_fiber'         , 13.47e9, ...
    'c22_r_fiber'       , 13.21e9, ...
    'c22_i_fiber'       , 0.49e9, ...
    'c44_fiber'         , 3.46e9, ...
    'c44_r_fiber'       , 3.29e9, ...
    'c44_i_fiber'       , 0.150e9, ...
    'c55_fiber'         , 5.20e9, ...
    'c55_r_fiber'       , 5.15e9, ...
    'c55_i_fiber'       , 0.30e9, ...
    'c66_fiber'         , 5.21e9, ...
    'c66_r_fiber'       , 5.11e9, ...
    'c66_i_fiber'       , 0.30e9, ...
    'E_resin'           , 3.7e9, ...
    'poisson_resin'     , 0.4, ...
    'G_resin'           , 1.32e9, ...
    'rho_fiber'         , 1528, ...  % density
    'rho_resin'         , 1270, ... % to be determined further
    'db_fiber'          , 0.1e-3, ... % complex velcity includes the info. of damping
    'db_resin'          , 0.15e-3, ...
    'h_fiber'           , 220e-6, ... % thinknesses
    'h_resin'           , 10e-6, ...
    'alpha_fiber_shear' , 0, ... % whatever, not used ...
    'alpha_resin_shear' , 0, ...
    'v_fiber_shear'     , 0, ...
    'v_resin_shear'     , 0, ...
    'elastic_matrix'    , elastic_constants_60per, ...
    'vf'                , 0.60, ...
    'rho_purefiber'     , 1800, ...
    'elastic_purefiber' , elastic_constants, ...
    'h_resin_ratio'     , 0); % set elastic_matrix as empty array [] to not to use it.

%%
elastic_matrix_ply      = Sample.elastic_purefiber;
% E11 longitudinal
elastic_matrix_ply(1)   = elastic_matrix_ply(1) .* Sample.vf + Sample.E_resin * (1 - Sample.vf);
% E22
elastic_matrix_ply(2:3) = Sample.E_resin ./ ( 1 - sqrt(Sample.vf) .* (1 - Sample.E_resin ./ elastic_matrix_ply(2:3)) );
% v12
elastic_matrix_ply(4:6) = elastic_matrix_ply(4:6) .* Sample.vf + Sample.poisson_resin * (1 - Sample.vf);
% G12
elastic_matrix_ply(7:9) = Sample.G_resin ./ ( 1 - sqrt(Sample.vf) .* (1 - Sample.G_resin ./ elastic_matrix_ply(7:9)) );
% E22 and E33 transverse
%             elastic_matrix_ply(2:3) =  ( Sample.vf ./ elastic_matrix_ply(2:3) + (1 - Sample.vf) ./ Sample.E_resin).^(-1);
%             elastic_matrix_ply(2:3) = elastic_matrix_ply(2:3) .* Sample.vf + Sample.E_resin * (1 - Sample.vf);
% possion's ratio v12 v13
elastic_matrix_ply(4:5) = elastic_matrix_ply(4:5) .* Sample.vf + Sample.poisson_resin * (1 - Sample.vf);
% G12 G13
elastic_matrix_ply(7:8) =  1 ./ ( Sample.vf ./  elastic_matrix_ply(7:8) + (1 - Sample.vf) ./ Sample.G_resin );

stiffness_matrix = fx_elastic2stiffness(elastic_matrix_ply);
Sample.c22_fiber = stiffness_matrix(2, 2);
Sample.c44_fiber = stiffness_matrix(4, 4);
Sample.c55_fiber = stiffness_matrix(5, 5);
Sample.c66_fiber = stiffness_matrix(6, 6);

% rule of mixture
rho_ply = Sample.vf .* Sample.rho_purefiber + (1 - Sample.vf) .* Sample.rho_resin;
v_fiber = sqrt(Sample.c22_fiber / rho_ply);

%% rotate the matrix

phi   = [0 pi/4 pi/2 pi*3/4];
theta = [0 0    0    0];
psi   = [0 0    0    0];

mat_paras = nan(length(psi), 22);

for i = 1:length(psi)
    stiffness_matrix_rotate = fx_StiffnessMatrixRotate2(stiffness_matrix, ...
        psi(i), theta(i), phi(i));
    
    %     stiffness_para = [...
    %         stiffness_matrix_rotate(1,1); stiffness_matrix_rotate(1:2, 2); ...
    %         stiffness_matrix_rotate(1:3, 3); stiffness_matrix_rotate(1:4, 4); ...
    %         stiffness_matrix_rotate(1:5, 5); stiffness_matrix_rotate(1:6, 6)]';
    stiffness_para = [...
        stiffness_matrix_rotate(1,1); stiffness_matrix_rotate(1:2, 2); ...
        stiffness_matrix_rotate(1:3, 3); stiffness_matrix_rotate(1:4, 4); ...
        stiffness_matrix_rotate(1:5, 5); stiffness_matrix_rotate(1:6, 6)]';
    mat_paras(i, :) = [stiffness_para rho_ply];
end

% save the material properties to a text file
writematrix(mat_paras, 'anisotroic_material_prop.txt','Delimiter','tab')

%%
% % Function to calculate the 3x3 rotation matrix
% function R = eulerToRotationMatrix(psi, theta, phi)
% Rz_psi = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
% Rx_theta = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% Rz_phi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
% 
% R = Rz_psi * Rx_theta * Rz_phi;
% end
% 
% % Function to calculate the 6x6 transformation matrix
% function T = rotationToTransformationMatrix(R)
% T = [R(1,1)^2, R(2,1)^2, R(3,1)^2, 2*R(2,1)*R(3,1), 2*R(3,1)*R(1,1), 2*R(1,1)*R(2,1);
%     R(1,2)^2, R(2,2)^2, R(3,2)^2, 2*R(2,2)*R(3,2), 2*R(3,2)*R(1,2), 2*R(1,2)*R(2,2);
%     R(1,3)^2, R(2,3)^2, R(3,3)^2, 2*R(2,3)*R(3,3), 2*R(3,3)*R(1,3), 2*R(1,3)*R(2,3);
%     R(1,2)*R(1,3), R(2,2)*R(2,3), R(3,2)*R(3,3), R(2,2)*R(3,3)+R(2,3)*R(3,2), R(3,2)*R(1,3)+R(3,3)*R(1,2), R(1,2)*R(2,3)+R(1,3)*R(2,2);
%     R(1,3)*R(1,1), R(2,3)*R(2,1), R(3,3)*R(3,1), R(2,3)*R(3,1)+R(2,1)*R(3,3), R(3,3)*R(1,1)+R(3,1)*R(1,3), R(1,3)*R(2,1)+R(1,1)*R(2,3);
%     R(1,1)*R(1,2), R(2,1)*R(2,2), R(3,1)*R(3,2), R(2,1)*R(3,2)+R(2,2)*R(3,1), R(3,1)*R(1,2)+R(3,2)*R(1,1), R(1,1)*R(2,2)+R(1,2)*R(2,1)];
% end
% 
% % Function to rotate the stiffness matrix
% function C_rotated = rotateStiffnessMatrix(C, psi, theta, phi)
% R = eulerToRotationMatrix(psi, theta, phi);
% T = rotationToTransformationMatrix(R);
% C_rotated = T' * C * T;
% end
% 



