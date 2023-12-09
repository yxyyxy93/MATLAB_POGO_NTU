% The code is intended for 
clear;
close all;
clc
addpath(genpath('tools'));
filename = 'E:\Backup\PhDdoc\Emerson\UIT20230208\EBSD\Polycrystals_POGOsimu\EBSD_data\EBSD_1521_V5_1.txt';
% filename = 'E:\Backup\PhDdoc\Emerson\UIT20230208\EBSD\Water2Polycrystals_POGOsimu\EBSD_data\EBSD_1521_V7_5.txt';
% filename = 'E:\Backup\PhDdoc\Emerson\UIT20230208\EBSD\Water2Polycrystals_POGOsimu\EBSD_data\EBSD_1532_V9.txt';
% filename = 'E:\Backup\PhDdoc\Emerson\UIT20230208\EBSD\Water2Polycrystals_POGOsimu\EBSD_data\EBSD_1521_V5_1_5um.txt';
% filename = 'E:\Backup\PhDdoc\Emerson\UIT20230208\EBSD\Water2Polycrystals_POGOsimu\EBSD_data\EBSD_1521_V7_5_5um.txt';
% filename = 'E:\Backup\PhDdoc\Emerson\UIT20230208\EBSD\Water2Polycrystals_POGOsimu\EBSD_data\EBSD_1532_V9_5um.txt';

% PogoFilename= 'EBSD_1521_V5_1_5um'; 
% PogoFilename= 'EBSD_1521_V7_5_5um'; 
% PogoFilename= 'EBSD_1532_V9_5um'; 
PogoFilename= 'EBSD_1521_V5_1_StrainMap'; 
% PogoFilename= 'EBSD_1521_V7_5_AttMap'; 
% % PogoFilename= 'EBSD_1532_V9_AttMap'; 

model.grain_Orientations = dlmread(filename)';

% X_sample = 1.8e-3;%unit: m.
% Y_sample = 1.8e-3;
X_mesh = 10e-6; % for 10MHz, near zone length is 2.43 mm, and 3.65 for 15 MHz
Y_mesh = 10e-6;
Z_mesh = 10e-6;
% ele_x = round(X_sample/X_mesh);
% ele_y = round(Y_sample/Y_mesh);


%% Stimulation signal
frequency = 15000;% unit: kHz
cycles    = 3;
timedelay = 0;
timestep  = 2e-10;
endtime   = 2e-6;
phase     = 0;
filename  = ['tb_',num2str(frequency),'kHz_',num2str(cycles),'cyc_',num2str(timestep),'_',num2str(endtime),'.dat'];
tb_signal=tbgeneration_sun(frequency, cycles, timedelay, timestep, endtime,phase,filename);
% tb_signal=continuousSin(frequency, timestep, endtime,phase);

%% Generating cubic and meshes
% nodes
X_cubic = 1.8e-3;%unit: m.
Y_cubic = 1.8e-3;
Z_cubic = Z_mesh;
X_nodePos = 0:X_mesh:X_cubic;
Y_nodePos = 0:Y_mesh:Y_cubic;
Z_nodePos = 0:Z_mesh:Z_cubic;
node_num = length(X_nodePos)*length(Y_nodePos)*length(Z_nodePos);

% for z = 1:length(Z_nodePos)
%     for y = 1:length(Y_nodePos)
%         for x = 1:length(X_nodePos)
%             model.nodePos(1,x+(y-1)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos)) = X_nodePos(x);
%             model.nodePos(2,x+(y-1)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos)) = Y_nodePos(y);
%             model.nodePos(3,x+(y-1)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos)) = Z_nodePos(z);
%         end
%     end
% end
% 
% % elements
% for z = 1:length(Z_nodePos)-1
%     for y = 1:length(Y_nodePos)-1
%         for x = 1:length(X_nodePos)-1
%             model.elNodes(1,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y-1)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos);
%             
%             model.elNodes(2,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y-1)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos)+1;
%             
%             model.elNodes(3,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos);
%             
%             model.elNodes(4,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y)*length(X_nodePos)+(z-1)*length(X_nodePos)*length(Y_nodePos)+1;
%             
%             model.elNodes(5,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y-1)*length(X_nodePos)+(z)*length(X_nodePos)*length(Y_nodePos);
%             
%             model.elNodes(6,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y-1)*length(X_nodePos)+(z)*length(X_nodePos)*length(Y_nodePos)+1;
%             
%             model.elNodes(7,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y)*length(X_nodePos)+(z)*length(X_nodePos)*length(Y_nodePos);
%             
%             model.elNodes(8,x+(y-1)*(length(X_nodePos)-1)+(z-1)*(length(X_nodePos)-1)*(length(Y_nodePos)-1))...
%                 = x+(y)*length(X_nodePos)+(z)*length(X_nodePos)*length(Y_nodePos)+1;
%         end
%     end
% end
% Lengths
lx = length(X_nodePos);
ly = length(Y_nodePos);
lz = length(Z_nodePos);
tic;
% 3D grids
[X, Y, Z] = ndgrid(X_nodePos, Y_nodePos, Z_nodePos);
% Node positions
model.nodePos = [X(:)'; Y(:)'; Z(:)'];
% Element nodes
[X, Y, Z] = ndgrid(1:(lx-1), 1:(ly-1), 1:(lz-1));
idx = @(x, y, z) x + (y-1)*lx + (z-1)*lx*ly;
model.elNodes = [reshape(idx(X(:), Y(:), Z(:)), 1, []); 
reshape(idx(X(:)+1, Y(:), Z(:)), 1, []); 
reshape(idx(X(:)+1, Y(:)+1, Z(:)), 1, []);
reshape(idx(X(:), Y(:)+1, Z(:)), 1, []);
reshape(idx(X(:), Y(:), Z(:)+1), 1, []); 
reshape(idx(X(:)+1, Y(:), Z(:)+1), 1, []);
reshape(idx(X(:)+1, Y(:)+1, Z(:)+1), 1, []);
reshape(idx(X(:), Y(:)+1, Z(:)+1), 1, [])];
toc;
model.elTypes{1}.name = 'C3D8';
model.elTypes{1}.paramsType = 0;
model.nDims = 3;
model.nDofPerNode = 3;
model.elTypeRefs = ones(length(model.elNodes(1,:)),1);

%% Put sample into a istropic metal for adding absorbing region
% model.matTypeRefs = ones(length(model.elNodes(1,:)),1);
% center = ones(2,length(model.elNodes(1,:)));
% for i=1:length(model.elNodes(1,:))
% 
%     center(1,i) = (model.nodePos(1,model.elNodes(1,i))+model.nodePos(1,model.elNodes(2,i))+...
%         model.nodePos(1,model.elNodes(3,i))+model.nodePos(1,model.elNodes(4,i)))/4;
%     
%     center(2,i) = (model.nodePos(2,model.elNodes(1,i))+model.nodePos(2,model.elNodes(2,i))+...
%         model.nodePos(2,model.elNodes(3,i))+model.nodePos(2,model.elNodes(4,i)))/4;
% end
% 
% X_lim_low = (X_cubic-X_sample)/2;
% X_lim_up = X_lim_low+X_sample;
% Y_lim_low = (Y_cubic-Y_sample)/2;
% Y_lim_up = Y_lim_low+Y_sample;
% 
% ele_index=find(...
%               (center(1,:)>=X_lim_low)&...
%               (center(1,:)<=X_lim_up)&...
%               (center(2,:)>=Y_lim_low)&...
%               (center(2,:)<=Y_lim_up));
% for i = 2:length(ele_index)+1
%     model.matTypeRefs(ele_index(i-1),1) = i;
% end
% 
% Z=reshape(model.matTypeRefs,[round(X_cubic/X_mesh),round(Y_cubic/Y_mesh)])';
% figure
% imshow(Z/max(max(Z)))
% set(gca, 'YDir', 'normal')


if ispc
    figure;
    node_index_temp=randsample(size(model.nodePos,2),1000);
    node_samples=model.nodePos(:,node_index_temp);
    if model.nDims==3%polt
        scatter3(node_samples(1,:),node_samples(2,:),node_samples(3,:));
    else
        scatter(node_samples(1,:),node_samples(2,:));
    end
    axis tight;axis equal;
end
%% Settings except for material
model.prec=8;% Precision
model.runName='Job';
model.nt=round(endtime/timestep);
model.dt=timestep;

%% Generator
model.shots{1, 1}.ntSig=length(tb_signal);
model.shots{1, 1}.dtSig=tb_signal(2,1)-tb_signal(1,1);
node_index_generator = [];
% Plane wave
Y_pos = Y_cubic;
Z_pos = Z_cubic;

X_lim_up  = X_cubic+X_mesh/2;
X_lim_low = 0-X_mesh/2;
Y_lim_up  = Y_pos+Y_mesh/2;
Y_lim_low = Y_pos-Y_mesh/2;
Z_lim_up  = Z_pos+Z_mesh/2;
Z_lim_low = 0-Z_mesh/2;
node_index_generator=find(...
                (model.nodePos(1,:)>=X_lim_low)&...
                (model.nodePos(1,:)<=X_lim_up)&...
                (model.nodePos(2,:)>=Y_lim_low)&...
                (model.nodePos(2,:)<=Y_lim_up)&...
                (model.nodePos(3,:)>=Z_lim_low)&...
                (model.nodePos(3,:)<=Z_lim_up)...
                );

hold on;scatter(model.nodePos(1,node_index_generator),model.nodePos(2,node_index_generator));
model.shots{1, 1}.sigs{1, 1}.sigType=1;% 0 - force, 1 - displacement
model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(node_index_generator),1)*2;
model.shots{1, 1}.sigs{1, 1}.nodeSpec=node_index_generator';
model.shots{1, 1}.sigs{1, 1}.sigAmps=ones(length(model.shots{1}.sigs{1}.dofSpec),1)*1e-13;
model.shots{1, 1}.sigs{1, 1}.sig=tb_signal(:,2);

%% Boundary
X_lim_up = 0;
X_lim_low = 0;
Y_lim_up = Y_cubic;
Y_lim_low = 0;
Z_lim_up = Z_cubic;
Z_lim_low = 0;
node_index_boundary = find(...
            (model.nodePos(1,:)>=X_lim_low)&...
            (model.nodePos(1,:)<=X_lim_up)&...
            (model.nodePos(2,:)>=Y_lim_low)&...
            (model.nodePos(2,:)<=Y_lim_up)&...
            (model.nodePos(3,:)>=Z_lim_low)&...
            (model.nodePos(3,:)<=Z_lim_up)...
            );
model.fixNodes = node_index_boundary;
model.fixDof=ones(length(node_index_boundary),1)*1;

X_lim_up = X_cubic;
X_lim_low = 0;
Y_lim_up = 0;
Y_lim_low = 0;
Z_lim_up = Z_cubic;
Z_lim_low = 0;
node_index_boundary = find(...
            (model.nodePos(1,:)>=X_lim_low)&...
            (model.nodePos(1,:)<=X_lim_up)&...
            (model.nodePos(2,:)>=Y_lim_low)&...
            (model.nodePos(2,:)<=Y_lim_up)&...
            (model.nodePos(3,:)>=Z_lim_low)&...
            (model.nodePos(3,:)<=Z_lim_up)...
            );
model.fixNodes = [model.fixNodes node_index_boundary];
model.fixDof=[model.fixDof; ones(length(node_index_boundary),1)*2];

X_lim_up = X_cubic;
X_lim_low = 0;
Y_lim_up = Y_cubic;
Y_lim_low = 0;
Z_lim_up = 0;
Z_lim_low = 0;
node_index_boundary = find(...
            (model.nodePos(1,:)>=X_lim_low)&...
            (model.nodePos(1,:)<=X_lim_up)&...
            (model.nodePos(2,:)>=Y_lim_low)&...
            (model.nodePos(2,:)<=Y_lim_up)&...
            (model.nodePos(3,:)>=Z_lim_low)&...
            (model.nodePos(3,:)<=Z_lim_up)...
            );
model.fixNodes = [model.fixNodes node_index_boundary];
model.fixDof=[model.fixDof; ones(length(node_index_boundary),1)*3];

hold on;scatter3(model.nodePos(1,model.fixNodes),model.nodePos(2,model.fixNodes),model.nodePos(3,model.fixNodes));


%% Receiver
% node_index_receiver=node_index_generator; % pulse-echo
node_index_receiver=[];

X_lim_up = X_cubic+X_mesh/2;
X_lim_low = 0-X_mesh/2;
Y_lim_up = Y_cubic+Y_mesh/2;
Y_lim_low = 0-Y_mesh/2;
Z_lim_up = Z_cubic+Z_mesh/2;
Z_lim_low = 0-Z_mesh/2;
node_index_receiver =  find(...
            (model.nodePos(1,:)>=X_lim_low)&...
            (model.nodePos(1,:)<=X_lim_up)&...
            (model.nodePos(2,:)>=Y_lim_low)&...
            (model.nodePos(2,:)<=Y_lim_up)&...
            (model.nodePos(3,:)>=Z_lim_low)&...
            (model.nodePos(3,:)<=Z_lim_up)...
            );

% hold on;scatter3(model.nodePos(1,node_index_receiver),model.nodePos(2,node_index_receiver),model.nodePos(3,node_index_receiver));
model.measSets{1, 1}.name='main';
model.measSets{1, 1}.isDofGroup=0;
model.measSets{1, 1}.measDof=repmat((1:model.nDims)',length(node_index_receiver),1);
model.measSets{1, 1}.measNodes=reshape(repmat(node_index_receiver,model.nDims,1),length(node_index_receiver)*model.nDims,1);
model.measFreq=1;
model.measStart=1;
model.fieldStoreIncs=round((1:1:30)/60*model.nt)';
%% Material settings
% For isotropic material
%-----------------------------------------------------%
% cWater = 1480; rhoWater = 1000; visc=0;
% model.matTypes{1,1}.paramsType = 5;%water media   
% model.matTypes{1,1}.paramValues= [cWater, rhoWater, visc];
% model.matTypes{1,1}.paramsType=0;%solid media (Steel)
% model.matTypes{1,1}.paramValues=[200e9,0.3,8000];
%-----------------------------------------------------%
% For anisotropic material
%-----------------------------------------------------%
% model.grain_Orientations=zeros(3,100);
% C11=243.1e9;C12=138.1e9;C44=121.9e9;%Fe
% % C11=234.6e9;C12=145.4e9;C44=126.2e9;%Inconel
% ES=[C11 C12 C11 C12 C12 C11 0 0 0 C44 0 0 0 0 C44 0 0 0 0 0 C44];%21 elastic constant for cubic FE
% SM_origin=[ES(1) ES(2) ES(4) ES(7) ES(11) ES(16);
%            ES(2) ES(3) ES(5) ES(8) ES(12) ES(17);
%            ES(4) ES(5) ES(6) ES(9) ES(13) ES(18);
%            ES(7) ES(8) ES(9) ES(10) ES(14) ES(19);
%            ES(11) ES(12) ES(13) ES(14) ES(15) ES(20);
%            ES(16) ES(17) ES(18) ES(19) ES(20) ES(21);];
C_11 = 2.046e11;
C_12 = 1.377e11;
C_13 = 1.377e11;
C_21 = C_12;
C_22 = C_11;
C_23 = C_12;
C_31 = C_13;
C_32 = C_23;
C_33 = C_11;
C_44 = 1.262e11;
C_55 = C_44;
C_66 = C_11;

SM_origin = [C_11, C_12, C_13, 0, 0, 0;
             C_21, C_22, C_23, 0, 0, 0;
             C_31, C_32, C_33, 0, 0, 0;
             0, 0, 0, C_44, 0, 0;
             0, 0, 0, 0, C_55, 0;
             0, 0, 0, 0, 0, C_66];

for i=1:length(model.grain_Orientations(1,:))
    SM_rotated(:,:)=StiffnessMatrixRotate2(SM_origin,-model.grain_Orientations(1,i)*pi/180,-model.grain_Orientations(2,i)*pi/180,-model.grain_Orientations(3,i)*pi/180);
    ES_rotated(i,:)=SM_rotated(:)';
end
ES_rotated(:,30)=[];ES_rotated(:,23:24)=[];ES_rotated(:,16:18)=[];ES_rotated(:,9:12)=[];ES_rotated(:,2:6)=[];

% For absorbing region case, matTypes starts from 2 (1 is for isotropic
% metal)
% Density=8000;
% for i=1:length(model.grain_Orientations(1,:))
%     model.matTypes{i+1,1}.paramsType=2;%2 for anistropic 
%     model.matTypes{i+1,1}.paramValues=[ES_rotated(i,:),Density];
% end

Density=8000;
for i=1:length(model.grain_Orientations(1,:))
    model.matTypes{i,1}.paramsType=2;%2 for anistropic
    model.matTypes{i,1}.paramValues=[ES_rotated(i,:),Density];
    model.matTypeRefs(i,1) = i;
end

[X,X1] = meshgrid(min(model.nodePos(1,:))+X_mesh/2:X_mesh:max(model.nodePos(1,:)));
[Y1,Y] = meshgrid(min(model.nodePos(2,:))+Y_mesh/2:Y_mesh:max(model.nodePos(2,:)));

k=1;
for i=1:length(Y)
    for j=1:length(X)
    Z(k) = model.matTypes{model.matTypeRefs((i-1)*length(X)+j),1}.paramValues(1);
    k=k+1;
    end
end
Z=reshape(Z,[length(X),length(Y)])';
figure
imshow(Z/max(max(Z)))
set(gca, 'YDir', 'normal')

% % Absorbing Regions
% % Only isotropic materials supported for SRM (stiffness reduction method)
% xLims = [0 0.9e-3 2.7e-3 3.6e-3];
% yLims = [0 0.9e-3 3.15e-3 3.6e-3];
% zLims = [];
% nAbsVals = 60;
% c0 = 5800 ;
% freq = 15e6;
% model = addAbsBound ( model, xLims , yLims , zLims , nAbsVals , [] , c0 , freq );
%% Save pogo-inp file
model = rmfield(model,'grain_Orientations');
savePogoInp(sprintf([PogoFilename,'_15MHz_C3D8.pogo-inp']),model, 1, 15);  % new version POGO
