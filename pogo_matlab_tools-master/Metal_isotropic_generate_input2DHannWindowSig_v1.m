% *NOTE:* althought the code has case for left right top bottom source
% position, but only TOP, pointSource is tested. For left/right/bottom DOF need to change.
% This code can generate isotropic and single crystal anisotropic POGO model.

% Generate pogo-input file for 2D single crystal (C) anisotropic or Isotropic model (E,nu,rho)
% Plane-Wave is considered
% unit: Hz,m, sec. (case sensitive)

%-------------------------------------------------------------------------
% UPDATES:
% update:- abhishek 11 April 2022, updated for pointsource and omnisource. But
% omni source is not correct because it is difficult to apply dispacement
% or stress equally in all direction.
% update:- abhishek 03 May 2022, added single crystal anisotropic material properties
% update:- abhishek 05 May 2022, Periodic BC (PBC)
% update:- abhishek 01 Sep 2022, displacement fluid
%-------------------------------------------------------------------------
% Set the following paramenters
% model Origin (0,0) is at Top left corner; Moving right x +(ve), moving down
% 'y' -(ve)
% (1) Excitation signal parameters (Source-time function);
% (2) Define model Dimensions and mesh size;
% (3) Choosee symmetry Boundary conditions if required,
% (4) Source and receivers
% (5) Material settings - isotropic (E,nu,rho)
%-------------------------------------------------------------------------
% NOTES:
% Regarding displacement fluid type.
% MeshElType = 'CPE4R';   signalType = 1; 

% INPUTS:
% PogoFilename      = File name
% symmBC            = symmetric boundary condition (BC)
% srcPos, recPos    =   source and receiver position
% AmpEndNodeHalf    =  corner nodes amplitude half for symmetric BC, if 
% generating plane wave
% absBound          = absorbing boundary
% PeriodicBC        = periodic BC
% materialType      = material type
% MeshElType        = mesh element type
% signalType        = excitation type.
% elesize           = mesh size
% Dim_x, Dim_y      = model dimension in x and y direction, following
% cartesian coordinate system.
% absRegionLenX, absRegionLenX = absorbing region dimension in x and y
% E, nu, rho        = Youngs modulus, Poissons ratio and density for
% isotropic media
% For single crystal anisotropic material properties
% epsilon, theta, si = Euler angles for anisotropic crystal rotation
% c11  to c66       = stiffness constants for anisotropic media
% freq, nCyc        = freq(Hz), number of cycles  		
% endtime           = totaime time for simulation  (endtime>2*Dim_y/cL)
% timestep          = timestep<0.8*elesize/cL; 
% measFreq          = measurement frequency
%%
close all;
clear; clc;
addpath(genpath('/pogo_matlab_tools-master/'));

%% Choose Parameter setting 
PogoFilename= 'JobAcousticIsotropic';   % name of input file to be saved
symmBC = 'no';                  % yes; no;
srcPos = 'top';                 % Source position for plane wave: top; bottom; left;right;pointSourceAtTop; planewave; pointSourceAtCenter; OmniSourceAtCenter;arbitrary
recPos = 'all';                 % Receivers: top; bottom; left;right; all; pointRecAtTop
AmpEndNodeHalf = 'no';          % 'yes'; 'no'; to keep the end-nodes half for Sym BC
absBound = 'no';                % yes; no;
PeriodicBC = 'PBC_L2R';         % 'PBC_T2B', 'PBC_L2R', 'PBC_all','noPBC'
materialType= 'isotropic';      % 'isotropic'; 'anisotropic'; 'acoustic'; 'dispFluid' 
MeshElType = 'CPE4R';           % 'CPE4'; CPE4R 'AC2D4';
signalType = 0;                 % 0, for force (recommended for solids), 1= displacement (displacement type fluid)

%% Define model parameters
elesize = 5e-6;

% model dimension
Dim_x = 15e-3;
Dim_y = 15e-3;
absRegionLenX = 0e-3;   % absorbing region length
absRegionLenY = 0e-3;

%% Setting model parameters
if strcmp(absBound,'no')
    absRegionLenX = 0e-3;
    absRegionLenY = 0e-3;
end

switch absBound
    case 'yes'
        absDim_x = absRegionLenX;
        absDim_y = absRegionLenY;
    case 'no'
        absDim_x = 0e-3;
        absDim_y = 0e-3;
end
        
totalDim_x = Dim_x + 2*absDim_x;    % absorbing boundary in both the sides of the model.
totalDim_y = Dim_y + 1*absDim_y;

% centre locations for grid, origin (0,0) at top left
% default the origin is at center of model if cx and cy =0
switch absBound
    case 'yes'
        cx = Dim_x/2;
        cy = -Dim_y/2;
    case 'no'
        cx = totalDim_x/2;
        cy = -totalDim_y/2;
end

dx=elesize;
dy=elesize;

nx=floor(totalDim_x/dx)+1;
ny=floor(totalDim_y/dx)+1;

%% Generate grid mesh
model = genGrid2D(nx,ny,dx,dx,cx,cy);
figure;
plot(model.nodePos(1,:),model.nodePos(2,:),'.');axis tight;axis equal;

%% Material properties: 
% For isotropic material
% G = 80e9; nu = E/2/G-1;
E=120e9;        %200e9;        %Youngs Modulus
nu=0.3;         %Poissons ratio
rho=4500;       %7800;       %Density
%sound speeds
cL = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)));
cSh = sqrt(E/(2*rho*(1+nu)));

% Acoustic: Water

% For single crystal anisotropic material properties
epsilon =0;
theta = 90;
psi= 90;

% Anisotropic: Hexagonal
c11=162e9;
c33=181e9;
c12=92e9;
c13=69e9;
c44=46.7e9;
SM_origin=[ c11,    c12,    c13,    0,      0,      0;
            c12,    c11,    c13,    0,      0,      0;
            c13,    c13,    c33,    0,      0,      0;
            0,      0,      0,      c44,    0,      0;
            0,      0,      0,      0,      c44,    0;
            0,      0,      0,      0,      0,      (c11-c12)/2 ];
Euler1_temp=deg2rad(epsilon);
Euler2_temp=deg2rad(theta);
Euler3_temp=deg2rad(psi);

%% Input (source-time function), output signal and filled data. 
freq = 10e6;                  % unit: Hz
nCyc = 3;                     % number of cycles  		
% Total_time=1e-5;      			
timestep = 2e-10;            % timestep<0.8*elesize/cL; 0.12*elesize/1500
endtime = 9e-6;             % endtime>2*Dim_y/cL
timePoints=round(endtime/timestep);

measFreq = 30; %5                      % measurement smapling frequency in output signal (number of increments between measurment being recorded)
gap = round(timePoints/80);           % do at around 80 times steps
if gap <1
    gap =1;
end
% fieldStoreIncs = 1:gap:timePoints;
fieldStoreIncs = 20;           % To view filed data. Field data storage increment: 50=more steps, 5=less steps
courant =0.3;
if timestep>(courant*elesize/cL)
    warning(['timestep must be less than ',num2str(elesize/cL*courant)])
elseif endtime<Dim_y/cL*2
    warning('Total time is not sufficient to capture all reflection')
end

% get wavelengths
lambda = cL/freq;
lambdaSh = cSh/freq;

if elesize > (lambda/20) 
    warning(['Element Size is not sufficient',num2str(0.8*elesize/cL)])
end
%% Symmetry BC, automatically define based on sources
switch symmBC
    case 'yes'
       if strcmp(srcPos,'top') | strcmp(srcPos,'bottom') 
            % Source at top, SymBdry at left and right
            model.fixNodes=find(...
                            (model.nodePos(1,:)==min(model.nodePos(1,:)))|...
                            (model.nodePos(1,:)==max(model.nodePos(1,:)))...
                            );
            model.fixDof=ones(size(model.fixNodes))*1;
        elseif strcmp(srcPos,'left') | strcmp(srcPos,'right')
            % % Source at left, SymBdry at top and bottom
            model.fixNodes=find(...
                (model.nodePos(2,:)==min(model.nodePos(2,:)))|...
                (model.nodePos(2,:)==max(model.nodePos(2,:)))...
                );
            model.fixDof=ones(size(model.fixNodes))*2;
       end
    case 'no'
        disp('No symmetric boundary condition');
end

%% Model settings except for material
model.prec=8;   % Precision
model.nDims=2;
model.nDofPerNode=2;
model.elTypes{1, 1}.name=MeshElType;
model.runName='Job';
%set dt and nt,want T sec long simulation with dt step
model.nt=timePoints;
model.dt=timestep;
% model.orientRefs;
% model.or;

%put a Hann signal on sig 1, frame 1
model = genPogoHannSignal(model,nCyc,freq,1,1);

%% Generator/Source:- Define the measurement sets.
switch srcPos
    case 'top'
        if strcmp(absBound,'yes')
            SrcNodeIndex = find(model.nodePos(2,:)==max(model.nodePos(2,:)) & ...
                model.nodePos(1,:)>=(Dim_x-Dim_x-absRegionLenX) & model.nodePos(1,:) <= (Dim_x+absRegionLenX));     % source at all nodes on top surface 
%                 model.nodePos(1,:)>=(Dim_x-Dim_x-absRegionLenX/2) & model.nodePos(1,:) <= (Dim_x+absRegionLenX/2)); % source partially at at absorbing region aslo 
%                 model.nodePos(1,:)>=(Dim_x-Dim_x) & model.nodePos(1,:) <=(Dim_x));     % source only at main model, not on absorbing region
        elseif strcmp(absBound,'no') 
            SrcNodeIndex=find(model.nodePos(2,:)==max(model.nodePos(2,:)));   % source on top at x=0,y=0
        end
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                               %set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*2;           %and which degrees of freedom
    case 'bottom'
        SrcNodeIndex=find(model.nodePos(2,:)==min(model.nodePos(2,:)));   % source on bottom 
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                                        %set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*2;           %and which degrees of freedom
    case 'left'
        SrcNodeIndex=find(model.nodePos(1,:)==min(model.nodePos(1,:)));   % source on left
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                                        %set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*1;           %and which degrees of freedom
    case 'right'
        SrcNodeIndex=find(model.nodePos(1,:)==max(model.nodePos(1,:)));   % source on right 
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                                        %set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*1;           %and which degrees of freedom
    case 'pointSourceAtTop'
%         SrcNodeIndex=find(model.nodePos(2,:)==(model.nodePos(2,:)));   
%         SrcNodeIndex=find((model.nodePos(2,:)==0e-3));      % for specific position
        SrcNodeIndex=find(model.nodePos(2,:)==max(model.nodePos(2,:)));                 % source on top at x=0,y=0
        SrcNodeIndex = SrcNodeIndex(ceil(numel(SrcNodeIndex)/2));                       %choose the center node at top position
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                               %set the type to 0, for force (recommended for solids), 1= displacement
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*2;           %and which degrees of freedom
    case 'pointSourceAtCenter'
        [~, idx_srcX]=min(abs(model.nodePos(1,:)-Dim_x/2)); %model.nodePos(1,idx_srcX)
        [~, idx_srcY]=min(abs(model.nodePos(2,:)+Dim_y/2)); %model.nodePos(2,idx_srcY)
        SrcNodeIndex=find(model.nodePos(1,:)==model.nodePos(1,idx_srcX) & ...
            model.nodePos(2,:)==model.nodePos(2,idx_srcY));
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                               % set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*2;           %apply load in direction 1, i.e. x
    case 'OmniSourceAtCenter'
        [~, idx_srcX]=min(abs(model.nodePos(1,:)-Dim_x/2)); %model.nodePos(1,idx_srcX)
        [~, idx_srcY]=min(abs(model.nodePos(2,:)+Dim_y/2)); %model.nodePos(2,idx_srcY)
        cneternode = find(model.nodePos(1,:)==model.nodePos(1,idx_srcX) & ...
            model.nodePos(2,:)==model.nodePos(2,idx_srcY));
        SrcNodeIndex = [cneternode];
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                                        %set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=(repmat(SrcNodeIndex,1,2))';             %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=[ones(length(SrcNodeIndex),1)*1; ones(length(SrcNodeIndex),1)*2];            %apply load in direction 1, i.e. x
        
    case 'arbitrary'
        [~, idx_srcX]=min(abs(model.nodePos(1,:)-Dim_x/3)); %model.nodePos(1,idx_srcX)
        [~, idx_srcY]=min(abs(model.nodePos(2,:)+Dim_y/3)); %model.nodePos(2,idx_srcY)
        SrcNodeIndex=find(model.nodePos(1,:)==model.nodePos(1,idx_srcX) & ...
            model.nodePos(2,:)==model.nodePos(2,idx_srcY));
        model.shots{1, 1}.sigs{1, 1}.sigType=signalType;                                        %set the type to 0, for force
        model.shots{1, 1}.sigs{1, 1}.isDofGroup=0;
        model.shots{1, 1}.sigs{1, 1}.nodeSpec=SrcNodeIndex';                           %say what nodes to apply signal to
        model.shots{1, 1}.sigs{1, 1}.dofSpec=ones(length(SrcNodeIndex),1)*2;           %apply load in direction 1, i.e. x
end

%% End node amplituude half or not for accurate symmetric BC.
switch AmpEndNodeHalf
    case 'no'
        if strcmp(srcPos,'OmniSourceAtCenter')
            model.shots{1, 1}.sigs{1, 1}.sigAmps=ones(length(model.shots{1}.sigs{1}.dofSpec),1)*1;    % this is to consider multiple dof
        else
            model.shots{1, 1}.sigs{1, 1}.sigAmps=ones(length(model.shots{1}.sigs{1}.dofSpec),1)*1;    % this is to consider multiple dof
        end
        
%         model.shots{1, 1}.sigs{1, 1}.sigAmps=ones(length(SrcNodeIndex),1)*1;           %set the signal amplitudes in each direction
        % model.shots{1, 1}.sigs{1, 1}.sig=tb_signal(:,2);                             % if toneburst signal generate using external function
        
    case 'yes'
        % To keep left and right node amplitude half for symmetric BC.
        model.shots{1, 1}.sigs{1, 1}.sigAmps=ones(length(SrcNodeIndex),1)*1;         %set the signal amplitudes in each direction
        model.shots{1, 1}.sigs{1, 1}.sigAmps([1,end])=0.5;                             %set the signal amplitudes 1/2 for end nodes in Sym BC        
end

%% Periodic BC
switch PeriodicBC
    case 'PBC_T2B'
        model.mpiDivNum = 1;
        model.mpi = cell(1,1);
        %top to bottom
        model.mpi{1}.outNodes = [(1:nx)+1*nx (1:nx)+(ny-2)*nx];
        model.mpi{1}.inNodes = [(1:nx)+(ny-1)*nx (1:nx)+0*nx];
    case 'PBC_L2R'
        model.mpiDivNum = 1;
        model.mpi = cell(1,1);
        %left to right
        model.mpi{1}.outNodes = [2 + (0:ny-1)*nx, nx-1+(0:ny-1)*nx];
        model.mpi{1}.inNodes = [nx+(0:ny-1)*nx, 1+(0:ny-1)*nx];
        %NB can do both of these together but need to set up for corner nodes to be
        %copied
    case 'PBC_all'
%         Corner nodes are not considered.
        model.mpiDivNum = 1;
        model.mpi = cell(1,1);
        
        % topOutNodes = [(1:nx)+1*nx (1:nx)+(ny-2)*nx];
        % topInNodes = [(1:nx)+(ny-1)*nx (1:nx)+0*nx];
        % not considering edge nodes
        topOutNodes = [(2:nx-1)+1*nx (2:nx-1)+(ny-2)*nx];
        topInNodes = [(2:nx-1)+(ny-1)*nx (2:nx-1)+0*nx];

        LeftOutNodes = [2 + (2:ny-3)*nx, nx-1+(2:ny-3)*nx];
        LeftInNodes = [nx+(2:ny-3)*nx, 1+(2:ny-3)*nx];
        
        model.mpi{1}.outNodes = [topOutNodes LeftOutNodes];
        model.mpi{1}.inNodes = [topInNodes LeftInNodes];
     
    case 'noPBC'
        disp('No Periodic BC applied');
end

%% Receivers 
model.measSets{1, 1}.name='main';
model.measSets{1, 1}.isDofGroup=0;

switch recPos
    case 'top'
        RecNodeIndex=find(model.nodePos(2,:)==max(model.nodePos(2,:)));   % source on top at x=0,y=0
%         RecNodeIndex = RecNodeIndex - length(RecNodeIndex);            %  if need second layer
        model.measSets{1, 1}.measDof= 2*ones(1,length(RecNodeIndex));       % 2*ones(1,no_el); 
        model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
    case 'bottom'
        RecNodeIndex=find(model.nodePos(2,:)==min(model.nodePos(2,:)));   % source on bottom 
        model.measSets{1, 1}.measDof= 2*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el); 
        model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
    case 'left'
        RecNodeIndex=find(model.nodePos(1,:)==min(model.nodePos(1,:)));   % source on left 
        model.measSets{1, 1}.measDof= 1*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el); 
        model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
    case 'right'
        RecNodeIndex=find(model.nodePos(1,:)==max(model.nodePos(1,:)));   % source on right 
        model.measSets{1, 1}.measDof= 1*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el); 
        model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
    case 'all'
        if strcmp(srcPos,'top') | strcmp(srcPos,'bottom') 
            if strcmp(absBound,'yes')
                RecNodeIndex=[1:length(model.nodePos(2,:))];
                model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
                model.measSets{1, 1}.measDof= 2*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el);
%                 find(model.nodePos(2,:)==max(model.nodePos(2,:)) & ...
%                 model.nodePos(1,:)>=0e-3 & model.nodePos(1,:) <= Dim_x);  
            elseif strcmp(absBound,'no') 
                RecNodeIndex=[1:length(model.nodePos(2,:))];
                model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
                model.measSets{1, 1}.measDof= 2*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el);
            end 
        elseif strcmp(srcPos,'left') | strcmp(srcPos,'right')
            RecNodeIndex=[1:length(model.nodePos(2,:))];
            model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
            model.measSets{1, 1}.measDof= 1*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el);
        elseif strcmp(srcPos,'pointSourceAtCenter')
            RecNodeIndex=[1:length(model.nodePos(2,:))];
            model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
            model.measSets{1, 1}.measDof= 2*ones(1,length(RecNodeIndex));  % 2*ones(1,no_el);
%             model.measSets{1, 1}.measDof= [1,2];  % 2*ones(1,no_el);
        elseif strcmp(srcPos,'OmniSourceAtCenter')
            RecNodeIndex=[1:length(model.nodePos(2,:))];
            model.measSets{1, 1}.measNodes=repelem(RecNodeIndex,2);  % nodeArray; % receive at these nodes
            model.measSets{1, 1}.measDof= (repmat([1,2],1,length(RecNodeIndex)));     %receive in x or y at node points         
        end
    case 'pointRecAtTop'
        RecNodeIndex = find(model.nodePos(2,:)==max(model.nodePos(2,:)));                 % source on top at x=0,y=0
        RecNodeIndex = RecNodeIndex(ceil(numel(RecNodeIndex)/2));
        model.measSets{1, 1}.measDof= 2*ones(1,length(RecNodeIndex));       % 2*ones(1,no_el); 
        model.measSets{1, 1}.measNodes=RecNodeIndex;  % nodeArray;
end

%% Measurement sets: define the measurement frequency and start time
model.measFreq=measFreq;
model.measStart=1;

%calculate time increments to store the field at
model.fieldStoreIncs=round((1:fieldStoreIncs)/fieldStoreIncs*model.nt)';     % to view field data
%% Material settings: 
% Type 0=Isotropic; 1=orthotropic;2=Anistropic;3=Engineering constants;
%4=Acoustic material; 5=displacement based fluid (experimental); 6=
%isotropic elastic and viscosity (experiemtnal)

if strcmp(materialType,'isotropic')
    % For isotropic material E, nu, rho, alpha (optional)
    model.matTypeRefs(1:length(model.elNodes),:)=1; %define the material references (set them all to material 1)
    model.matTypes{1,1}.paramsType=0;               %Type 0=Isotropic; 1=orthotropic;2=Anistropic
    model.matTypes{1,1}.paramValues= [E, nu, rho];
elseif strcmp(materialType,'anisotropic')
    % For an-isotropic material
    SM_rotated_all(:,:,1,1)=StiffnessMatrixRotate2(SM_origin,-Euler3_temp,-Euler2_temp,-Euler1_temp);
    model.matTypeRefs(1:length(model.elNodes),:)=1; %define the material references (set them all to material 1)
    model.matTypes{1,1}.paramsType=2;               %Type 0=Isotropic; 1=orthotropic;2=Anistropic
    C=SM_rotated_all;
    % Assigning material type
    model.matTypes{1,1}.paramValues=[C(1,1);C(1,2);C(2,2);C(1,3);C(2,3);C(3,3);C(1,6);C(2,6);...
                                        C(3,6);C(6,6);C(1,5);C(2,5);C(3,5);C(6,5);C(5,5);C(1,4);...
                                        C(2,4);C(3,4);C(6,4);C(5,4);C(4,4);...
                                        rho; 0];
elseif strcmp(materialType,'acoustic')
%   Material properties: c,rho,alpha (optional)
%  Material, name=water; *Acoustic Medium 2.2e+9, *Density 1000
    cWater = 1500; rhoWater = 1000; visc=0;
    model.matTypeRefs(1:length(model.elNodes),:)=1; %define the material references (set them all to material 1)
    model.matTypes{1,1}.paramsType = 4;               
    model.matTypes{1,1}.paramValues= [cWater, rhoWater, visc];


elseif strcmp(materialType,'dispFluid') % displacment based fluid
%   Material properties: c,rho,mu (viscosity), alpha (optional)
	cWater = 1500; rhoWater = 1000; visc=0;
    model.matTypeRefs(1:length(model.elNodes),:)=1; %define the material references (set them all to material 1)
    model.matTypes{1,1}.paramsType=5;               
    model.matTypes{1,1}.paramValues= [cWater, rhoWater, visc];    
end
%% Save pogo-inp file
% savePogoInp( [PogoFilename,'.pogo-inp'], model, 1.07 );
savePogoInp(sprintf([PogoFilename,'.pogo-inp']),model, 1, 15);  % new version POGO
%% Absorbing Region
% x_boundary_ROI=[0e-3,Dim_x];
% 
% switch absBound
%     case 'yes'
%         % Create absorbing layer 0e-3,Dim_x
%         xLims=[min(model.nodePos(1,:)),x_boundary_ROI(1),x_boundary_ROI(2),max(model.nodePos(1,:))];
%         yLims = []; 
%         zLims = [];
% %         yLims=[min(model.nodePos(2,:)),y_boundary_ROI(1),y_boundary_ROI(2),max(model.nodePos(2,:))];
%         nAbsVals=60;
% %         clearvars -except model xLims yLims zLims nAbsVals foldername absBound
% %          model = addAbsBound(model,xLims,yLims,zLims,[],[],cL,freq);
%          modelabs = addAbsBound(model,xLims,yLims,zLims, nAbsVals, [], 6131, 10e6, 0 );
%         
% %         Writing new pogo-inp file
%         savePogoInp( [PogoFilename,'_Absorb','.pogo-inp'], modelabs, 1.07 );
% 
%      case 'no'
%          disp('No absorbing region');
% end

%% visualisation
px = model.nodePos(1,:); 
py = model.nodePos(2,:);
figure;
plot(px,py,'.b','Markersize',3,'linewidth',1);hold on;
plot(px(RecNodeIndex),py(RecNodeIndex),'*k','Markersize',8,'linewidth',0.5);
plot(px(SrcNodeIndex),py(SrcNodeIndex),'og','Markersize',4,'linewidth',1,'MarkerFaceColor',[0.2 1 0.2]);
if strcmp(symmBC,'yes')
    plot(px(model.fixNodes),py(model.fixNodes),'.g','Markersize',8,'linewidth',1); 
    axis tight;axis equal;legend('model','receiver','source','symmetric BC');hold off;
elseif strcmp(PeriodicBC,'PBC_T2B') |  strcmp(PeriodicBC,'PBC_all') | strcmp(PeriodicBC,'PBC_L2R')
    hold on
    plot(px(model.mpi{1}.outNodes),py(model.mpi{1}.outNodes),'or')
    plot(px(model.mpi{1}.inNodes),py(model.mpi{1}.inNodes),'ok')
    axis tight;axis equal;legend('model','receiver','source','Periodic BC');hold off;
else
    axis tight;axis equal;legend('model','receiver','source');hold off;
end

% Add data labels to scatters
% figure;
% plot(px,py,'.b','Markersize',3,'linewidth',1);
% b=num2str(px); c=cellstr(b);
% text(px,py, b)

% x=1:10;
% y=1:10;
% scatter(x,y);
% a = [1:10]'; b=num2str(a); c=cellstr(b);
% 
% dx = 0.1; dy =0.1;
% text(x+dx, y+dy, c, 'Fontsize', 8);
%% Running Solver

% clear;clc; close all;
% filename = 'structured_2d';
% system(sprintf(['pogoBlock ',filename,'.pogo-inp']))
% system(sprintf(['pogoSolve ',filename,'.pogo-inp']))
% clear; close all;%clc;

%% Field Visualisation
% f = loadPogoField('structured_2d.pogo-field');
% % viewPogoField( f);
% frame_no = 88;
% npix=200; %default =100
% viewPogoField( f, npix, frame_no);colorbar;colormap jet;%caxis([-3e-11 3e-11]);
