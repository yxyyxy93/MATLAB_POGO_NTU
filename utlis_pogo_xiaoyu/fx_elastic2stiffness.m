function stiffness_matrix = fx_elastic2stiffness(elastic_constants)
% input: elastic_constants = [E1, E2, E3, v12, v13, v23, G12, G13, G23]
% unit: GPa

%% the output is the inverted compliance matrix !!!!!!
E1 = elastic_constants(1);
E2 = elastic_constants(2);
E3 = elastic_constants(3);

v12 = elastic_constants(4);
v13 = elastic_constants(5);
v23 = elastic_constants(6);

G12 = elastic_constants(7);
G13 = elastic_constants(8);
G23 = elastic_constants(9);

% transver isotropic
v21 = v12 * E2 / E1;
v32 = v23 * E3 / E2;
v31 = v13 * E3 / E1;

elastic_matrix_lt = [
    1/E1       -v21/E2 -v31/E3;
    -v12/E1 1/E2       -v32/E3;
    -v13/E1 -v23/E2 1/E3;
    ];

elastic_matrix_rb = [
    1 /  G23  0                  0;
    0                 1 / G13   0;
    0                  0                 1 / G12];

elastic_matrix = [
elastic_matrix_lt zeros(3);
zeros(3)               elastic_matrix_rb];

stiffness_matrix = inv(elastic_matrix);


end

