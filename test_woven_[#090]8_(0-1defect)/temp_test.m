currentFileFolder = "./base_model_shiftseed_1";
model_path        = "./base_model_shiftseed_4.mat";

% fx_Cscanread(currentFileFolder)
%
nx = 20;
ny = 20;
dx = 0.2e-3;
dy = 0.2e-3;

fx_immersion_Cscan(currentFileFolder, model_path, nx, ny, dx, dy);