clear;
close all;
clc
fclose all;

PogoFilename = 'immersionUT_CaseStudy_yxy'; 

% Stimulation signal
frequency = 5e6;% unit: Hz
cycles    = 3;
timedelay = 2e-7;
timestep  = 2e-10;
endtime   = 1e-6;
phase     = -pi/2;
filename  = ['tb_',num2str(frequency),'MHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal = tbgeneration_yxy(frequency, cycles, timedelay, timestep, endtime,phase,filename);

%% Generating cubic and meshes
% nodes
X_cubic   = 10e-3;  %unit: m.
Y_cubic   = 30e-3;
X_mesh    = 1.0e-5; % 20 elements per wavelength is minimum
Y_mesh    = 1.0e-5;
X_nodePos = 0:X_mesh:X_cubic;
Y_nodePos = 0:Y_mesh:Y_cubic;
node_num  = length(X_nodePos) * length(Y_nodePos);

for y = 1:length(Y_nodePos)
    for x = 1:length(X_nodePos)
        model.nodePos(1,x+(y-1) * length(X_nodePos)) = X_nodePos(x);
        model.nodePos(2,x+(y-1) * length(X_nodePos)) = Y_nodePos(y);
    end
end

% elements
for y = 1:length(Y_nodePos)-1
    for x = 1:length(X_nodePos)-1
        model.elNodes(1,x+(y-1)*(length(X_nodePos)-1))...
            = x+(y-1)*length(X_nodePos);
        
        model.elNodes(2,x+(y-1)*(length(X_nodePos)-1))...
            = x+(y-1)*length(X_nodePos)+1;
        
        model.elNodes(4,x+(y-1)*(length(X_nodePos)-1))...
            = x+(y)*length(X_nodePos);
        
        model.elNodes(3,x+(y-1)*(length(X_nodePos)-1))...
            = x+(y)*length(X_nodePos)+1;
    end
end

model.elTypes{1}.name       = 'CPE4R';
model.elTypes{1}.paramsType = 0;
model.nDims                 = 2; % default 2D 
model.nDofPerNode           = 2;
model.elTypeRefs            = ones(length(model.elNodes(1,:)),1);

% model.grain_Orientations(:,1) = [0 0 0]'; %<100>
% model.grain_Orientations(:,2) = [0 45 0]'; %<110>
% model.grain_Orientations(:,1) = [0 0 0]'; %<100>

% model.grain_Orientations(:,2) = [45 0 0]'; %<110>
% lamellar_width1 = 5e-3;%<110>
% lamellar_width2 = 5e-3;%<100>

model.matTypeRefs = ones(length(model.elNodes(1,:)), 1); %water medium

%% Select elements that in the range of solid parts and change their matTypeRefs
% to 2, which refers to the elastic media.

% Define the CFRP region
X_rect_start = 0;  % Starting X position of the rectangle
Y_rect_start = 20e-3; % Starting Y position of the rectangle
X_rect_end   = X_cubic; % Ending X position of the rectangle
Y_rect_end   = 20e-3+24*230e-6; % Ending Y position of the rectangle
model = fx_assign_material_to_rectangle(model, ...
    X_rect_start, Y_rect_start, X_rect_end, Y_rect_end, 2);

% Define interplies
X_rect_start = 0 * ones(1, 25); % Starting X position of the rectangle
Y_rect_start = 20e-3-10e-6: 230e-6: 20e-3+24*230e-6; % Starting Y position of the rectangle
X_rect_end   = X_cubic * ones(1, 25); % Ending X position of the rectangle
Y_rect_end   = 20e-3: 230e-6: 20e-3+24*230e-6; % Ending Y position of the rectangle
model = fx_assign_material_to_rectangle(model, ...
    X_rect_start, Y_rect_start, X_rect_end, Y_rect_end, 3);

pic = reshape(model.matTypeRefs, round(X_cubic/X_mesh), round(Y_cubic/Y_mesh));
figure;
imagesc(pic);

% if ispc
%     figure;
%     node_index_temp=randsample(size(model.nodePos,2),1000);
%     node_samples=model.nodePos(:,node_index_temp);
%     if model.nDims==3%polt
%         scatter3(node_samples(1,:),node_samples(2,:),node_samples(3,:));
%     else
%         scatter(node_samples(1,:),node_samples(2,:));
%     end
%     axis tight;axis equal;
% end

%% Settings except for material
model.prec    = 8;      % Precision
model.runName = 'Job';
model.nt      = 8e4;
model.dt      = 0.5e-9;

%% GenerateMicrograins
% if Dims == 2
% xx = [0.0250 0.0270 0.0265 0.0260 0.0255 0.0252 0.0250 0.0245 0.0240 0.0235 0.0230 0.0235 0.0240 0.0250];
% yy = [0.0180 0.0185 0.0150 0.0130 0.0090 0.0070 0.0065 0.0060 0.0055 0.0052 0.0050 0.0090 0.0130 0.0180];
% mean_EulerAngle = [0 0 0]';
% sigma_EulerAngle = [1 1 1]';
% num_of_added_grain = 50;
% alpha=1;
% beta=0;
% gamma=0;
% model = GenerateMicrograins( model, xx, yy, mean_EulerAngle, sigma_EulerAngle, num_of_added_grain, alpha, beta, gamma )
% end

%% Material settings
% For isotropi49c material
%-----------------------------------------------------%
cWater = 1500; rhoWater = 1000; visc = 0;

model.matTypes{1,1}.paramsType  = 5; % water media
model.matTypes{1,1}.paramValues = [cWater, rhoWater, visc];

model.matTypes{2,1}.paramsType  = 0; % CFRP
% model.matTypes{2,1}.paramValues = [200e9,0.3,7800];
model.matTypes{2,1}.paramValues = [13.47e9, 0, 1588]; % virtual parameter!

model.matTypes{3,1}.paramsType  = 0; % resin-rich interply
model.matTypes{3,1}.paramValues = [3.7e9, 0.4, 1270];


%-----------------------------------------------------%
% For anisotropic material
%-----------------------------------------------------%
model.grain_Orientations=zeros(3,100);

C11=243.1e9;C12=138.1e9;C44=121.9e9; %Fe
% C11=234.6e9;C12=145.4e9;C44=126.2e9; %Inconel

ES=[C11 C12 C11 C12 C12 C11 0 0 0 C44 0 0 0 0 C44 0 0 0 0 0 C44]; % 21 elastic constant for cubic FE

SM_origin=[ES(1) ES(2) ES(4) ES(7) ES(11) ES(16);
           ES(2) ES(3) ES(5) ES(8) ES(12) ES(17);
           ES(4) ES(5) ES(6) ES(9) ES(13) ES(18);
           ES(7) ES(8) ES(9) ES(10) ES(14) ES(19);
           ES(11) ES(12) ES(13) ES(14) ES(15) ES(20);
           ES(16) ES(17) ES(18) ES(19) ES(20) ES(21);];

% for i=1:length(model.grain_Orientations(1,:))
%     SM_rotated(:,:)=StiffnessMatrixRotate2(SM_origin,-model.grain_Orientations(1,i)*pi/180,-model.grain_Orientations(2,i)*pi/180,-model.grain_Orientations(3,i)*pi/180);
%     ES_rotated(i,:)=SM_rotated(:)';
% end
% ES_rotated(:,30)=[];ES_rotated(:,23:24)=[];ES_rotated(:,16:18)=[];ES_rotated(:,9:12)=[];ES_rotated(:,2:6)=[];

% matTypeRefs=model.matTypeRefs;
% model = rmfield(model,'matTypeRefs');
% model.matTypeRefs(1:length(model.elNodes),:)=matTypeRefs;

% Poissonratio=zeros(length(model.grain_Orientations),1);
% Poissonratio(:)=0.3;
% Density=zeros(length(model.grain_Orientations),1);
% Density(:)=7800;
% Density(:)=8260;
% for i=1:length(model.grain_Orientations(1,:))
% model.matTypes{i,1}.paramsType=2;%2 for anistropic
% model.matTypes{i,1}.paramValues=[ES_rotated(i,:),Density(i)];
% end

% figure
% X =
% elsize=Element_size;
Element_size = X_mesh;

[X, X1] = meshgrid(min(model.nodePos(1,:)) + Element_size/2: ...
    Element_size: max(model.nodePos(1,:)));
[Y1, Y] = meshgrid(min(model.nodePos(2,:)) + Element_size/2: ...
    Element_size:max(model.nodePos(2,:)));

% X = X(1:100,1:100);
% Y = Y(1:100,1:100);

for i=1:length(Y)
    for j=1:length(X)
        Z(i,j) = model.matTypes{model.matTypeRefs((i-1)*length(X)+j),1}.paramValues(3);
    end
end

[m, n]=size(Z);
X(m+1:end,:)=[];

for i=m+1:n
    Y(:,i)=Y(:,1);
end

mesh(X(1,:),Y(:,1)',Z);

colorbar;
xlabel('\itx\rm (mm)');
ylabel('\ity\rm (mm)');
zlabel('C22\rm(Pa)');
% axis([0 15e-3 0 15e-3])

% figure
% contour(X,Y,Z,100,'Linewidth',0.5);
% colorbar;
% % axis([0 15e-3 0 15e-3])
% xlabel('\itx\rm (mm)');
% ylabel('\ity\rm (mm)');

% Absorbing Regions
% Only isotropic materials supported for SRM (stiffness reduction method)
xLims    = [];
yLims    = [Y_cubic Y_cubic-2e-3 100 102];
zLims    = [];
nAbsVals = 60;
c0       = 1500 ;
freq     = frequency;
model    = addAbsBound(model, xLims , yLims , zLims , [], [], c0, freq);

%% Generator
model.shots{1, 1}.ntSig = length(tb_signal);
model.shots{1, 1}.dtSig = tb_signal(2,1)-tb_signal(1,1);
node_index_generator    = [];

% % Single transducer
% for i        = 16:18
%     Element_size = 10e-6;
%     X_pos        = 0;
%     Y_pos        = 0.002075+(i-1)*0.31e-3;
%     % X_lim_up   = max(model.nodePos(1,:));
%     % X_lim_low  = min(model.nodePos(1,:));
%     X_lim_up     = X_pos+Element_size/2;
%     X_lim_low    = X_pos-Element_size/2;
%     Y_lim_up     = Y_pos+Element_size/2;
%     Y_lim_low    = Y_pos-Element_size/2;
%     node_index_generator2=find(...
%         (model.nodePos(1,:)>=X_lim_low)&...
%         (model.nodePos(1,:)<=X_lim_up)&...
%         (model.nodePos(2,:)>=Y_lim_low)&...
%         (model.nodePos(2,:)<=Y_lim_up));
%     node_index_generator = [node_index_generator node_index_generator2];
% end

% Plane wave
X_lim_up  = 8e-3+Element_size/2;
X_lim_low = 2e-3-Element_size/2;
Y_lim_up  = 0+Element_size/2;
Y_lim_low = 0-Element_size/2;

node_index_generator = find( ...
    (model.nodePos(1,:) >= X_lim_low) & ...
    (model.nodePos(1,:) <= X_lim_up) & ...
    (model.nodePos(2,:) >= Y_lim_low) & ...
    (model.nodePos(2,:) <= Y_lim_up));

% Calculate local spacing
% [x_coor_sorted,index_sorted]=sort(model.nodePos(1,node_index_generator));
% LocalSpacing=(x_coor_sorted(3:end)-x_coor_sorted(1:end-2))/2;
% %             [y_coor_sorted,index_sorted]=sort(model.nodePos(2,node_index_generator));
% %             LocalSpacing=(y_coor_sorted(3:end)-y_coor_sorted(1:end-2))/2;
% LocalSpacing=5e-5;
% node_index_generator=node_index_generator(index_sorted(2:end-1));

hold on;
scatter(model.nodePos(1,node_index_generator), model.nodePos(2,node_index_generator));

fd     = 25.4e-3; % m
% diameter = 6.35e-3; 
center = [(X_lim_up+X_lim_low)/2 (Y_lim_low+Y_lim_up)/2];
focused_waves = fx_focused_wave(fd, center, timestep, ...
    model.nodePos(1:2, node_index_generator), cWater, tb_signal(:,2));

% figure;
% plot(focused_waves(1, :)');
% hold on;
% plot(focused_waves(301, :)');

for i = 1: size(focused_waves, 1)
    model.shots{1, 1}.sigs{i, 1}.sigType    = 0; % 0 - force, 1 - displacement
    model.shots{1, 1}.sigs{i, 1}.isDofGroup = 0;
    model.shots{1, 1}.sigs{i, 1}.dofSpec    = 2;
    model.shots{1, 1}.sigs{i, 1}.nodeSpec   = node_index_generator(i)';
    model.shots{1, 1}.sigs{i, 1}.sigAmps    = ones(length(model.shots{1}.sigs{i}.dofSpec),1)*1e-13;
    model.shots{1, 1}.sigs{i, 1}.sig        = focused_waves(i, :)';
end

%% Boundary
X_lim_up  = 0;
X_lim_low = 0;
Y_lim_up  = Y_cubic;
Y_lim_low = 0;
node_index_generator2=find(...
    (model.nodePos(1,:)>=X_lim_low)&...
    (model.nodePos(1,:)<=X_lim_up)&...
    (model.nodePos(2,:)>=Y_lim_low)&...
    (model.nodePos(2,:)<=Y_lim_up));
model.fixNodes = node_index_generator2;

X_lim_up  = X_cubic;
X_lim_low = X_cubic;
Y_lim_up  = Y_cubic;
Y_lim_low = 0;
node_index_generator2=find(...
    (model.nodePos(1,:)>=X_lim_low)&...
    (model.nodePos(1,:)<=X_lim_up)&...
    (model.nodePos(2,:)>=Y_lim_low)&...
    (model.nodePos(2,:)<=Y_lim_up));

model.fixNodes = [model.fixNodes node_index_generator2];

hold on;
scatter(model.nodePos(1,model.fixNodes),model.nodePos(2,model.fixNodes));
model.fixDof=ones(length(model.fixNodes),1)*1;

%% Receiver
node_index_receiver = node_index_generator;
% node_index_receiver=[];
% X_lim = linspace(16.4e-3,54.2e-3,64);
% for i = 1:length(X_lim)
% Element_size = 50e-6;
% X_pos = X_lim(i);
% Y_pos = 20e-3;
% X_lim_up = X_pos+Element_size/2;
% X_lim_low = X_pos-Element_size/2;
% Y_lim_up = Y_pos+Element_size/2;
% Y_lim_low = Y_pos-Element_size/2;
% node_index_receiver2=find(...
%                 (model.nodePos(1,:)>=X_lim_low)&...
%                 (model.nodePos(1,:)<=X_lim_up)&...
%                 (model.nodePos(2,:)>=Y_lim_low)&...
%                 (model.nodePos(2,:)<=Y_lim_up));
% node_index_receiver=[node_index_receiver node_index_receiver2];
% end
% X_lim_up = X_cubic/2+10e-6*3;
% X_lim_low = X_cubic/2-10e-6*3;
% Y_lim_up = 0;
% Y_lim_low = 0;
% Z_lim_up = Z_cubic/2 + 25e-6*3;
% Z_lim_low = Z_cubic/2 - 25e-6*3;
% node_index_receiver = find(...
%             (model.nodePos(1,:)>=X_lim_low)&...
%             (model.nodePos(1,:)<=X_lim_up)&...
%             (model.nodePos(2,:)>=Y_lim_low)&...
%             (model.nodePos(2,:)<=Y_lim_up)&...
%             (model.nodePos(3,:)>=Z_lim_low)&...
%             (model.nodePos(3,:)<=Z_lim_up)...
%             );
%
% X_lim_up = X_cubic/2+10e-6*3;
% X_lim_low = X_cubic/2-10e-6*3;
% Y_lim_up = Y_cubic;
% Y_lim_low = Y_cubic;
% Z_lim_up = Z_cubic/2+25e-6*3;
% Z_lim_low = Z_cubic/2-25e-6*3;
% node_index_receiver2 = find(...
%             (model.nodePos(1,:)>=X_lim_low)&...
%             (model.nodePos(1,:)<=X_lim_up)&...
%             (model.nodePos(2,:)>=Y_lim_low)&...
%             (model.nodePos(2,:)<=Y_lim_up)&...
%             (model.nodePos(3,:)>=Z_lim_low)&...
%             (model.nodePos(3,:)<=Z_lim_up)...
%             );
%
% node_index_receiver=[node_index_receiver node_index_receiver2];
hold on;
scatter(model.nodePos(1,node_index_receiver), model.nodePos(2,node_index_receiver));

model.measSets{1, 1}.name       = 'main';
model.measSets{1, 1}.isDofGroup = 0;
model.measSets{1, 1}.measDof    = repmat((1:model.nDims)',length(node_index_receiver),1);
model.measSets{1, 1}.measNodes  = reshape(repmat(node_index_receiver,model.nDims,1),length(node_index_receiver)*model.nDims,1);
model.measFreq       = 1;
model.measStart      = 1;
model.fieldStoreIncs = round((1:2:60)/60*model.nt)';

%% Save pogo-inp file
% model = rmfield(model,'grain_Orientations');
savePogoInp(sprintf([PogoFilename,'.pogo-inp']), model, 1, 15);  % new version POGO

disp(".pogo-inp saved");

% pic = reshape(model.matTypeRefs(:,1),X_cubic/X_mesh,Y_cubic/Y_mesh)-1;
% figure;
% imshow(pic');

close all;
