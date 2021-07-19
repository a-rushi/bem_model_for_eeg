clear all;
close all;
clc;

%% Important References
% http://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_dipole_fit/
% http://www.fieldtriptoolbox.org/reference/ft_prepare_leadfield/
% http://www.fieldtriptoolbox.org/reference/ft_prepare_sourcemodel/

%% Add path of fieldtrip toolbox
addpath('/Users/aarushigupta/Desktop/BEM_Amita/fieldtrip-20190224/fieldtrip-20190224')
ft_defaults  % initiate fieldtrip

%% Reading MRI file - anatomical data
% mri = ft_read_mri('standard_mri.mat');  % E:\EEG_DATA\fieldtrip-20190224\fieldtrip-20190224\template\headmodel\standard_mri.mat 
% cfg        = [];    % To plot MRI define configuration
% cfg.method = 'ortho'; % cfg.method could be 'slice', 'ortho', 'glassbrain'
% cfg.colorbar = 'no';
% ft_sourceplot(cfg, mri);  
% ft_determine_coordsys(mri, 'interactive', 'no');

%% Segmentation of the anatomical information into different tissue types
% cfg        = [];    % cfg: configuration of the function that was used to create vol
% cfg.output       = {'brain','skull','scalp'};
% segmentedmri  = ft_volumesegment(cfg, mri);
% % save segmentedmri segmentedmri

%% Define number of vertices of each segmented tissue to prepare triangulation mesh 
% load segmentedmri
% cfg        = [];
% cfg.tissue = {'brain','skull','scalp'};
% cfg.numvertices = [1500 1500 1500];
% bnd = ft_prepare_mesh(cfg,segmentedmri); % The bnd field contains information about the mesh
% % save bnd bnd

%% Prepare Headmodel - bemcp method i.e Boundary element method
% load bnd
% cfg = [];
% cfg.method = 'bemcp'; % You can also specify 'openmeeg', 'dipoli', or another method.
% cfg.conductivity = [0.3300 0.0041 0.3300];
% vol = ft_prepare_headmodel(cfg,bnd);
% % save vol vol

%% Visualization plot mesh
% load vol
% figure
% ft_plot_mesh(vol.bnd(1), 'facecolor',[0 0 1], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
% hold on;
% ft_plot_mesh(vol.bnd(2),'facecolor',[0 1 1],'edgecolor','none','facealpha',0.4);
% hold on;
% ft_plot_mesh(vol.bnd(3),'edgecolor','none','facealpha', 0.3,'facecolor',[0.4 0.6 0.4]);

%% Electrode placement over the scalp surface
% load vol
% elec = ft_read_sens('standard_1020.elc');  % for electrodes - 97 electrodes - default - 'standard_1020.elc'
% % save elec elec
% figure;
% ft_plot_mesh(vol.bnd(3), 'edgecolor','none','facealpha',0.3,'facecolor',[0.4 0.6 0.4]);
% hold on;
% ft_plot_sens(elec);     % optional properties - style,marker

%% Place source in prepared BEM head model, this step is called Source modeling
load vol
load elec
cfg = [];
cfg.grid.xgrid      = [5];    % source x position
cfg.grid.ygrid      = [-5];   % source y position
cfg.grid.zgrid      = [5];    % source z position
zplane = cfg.grid.zgrid;  
[sourcemodel] = ft_prepare_sourcemodel(cfg);
dip_pos = sourcemodel.pos; 
for i =1:length(dip_pos(:,1))
dip_posm(i,:) = dip_pos(i,:)/norm(dip_pos(i,:));  % radial dipole orientation
end
figure
ft_plot_mesh(vol.bnd(1),'edgecolor','none','facealpha',0.5,'facecolor',[0.6 0.6 0.8]);
hold on
ft_plot_dipole(dip_pos,dip_posm','thickness',3,'length',20,'diameter',10);

%% Computation of simulated EEG recording using the Fieldtrip function
cfg      = [];
cfg.headmodel = vol;          % see above
cfg.elec = elec;              % see above
cfg.dip.pos = dip_pos;
cfg.dip.mom = reshape(dip_posm',1,3*length(dip_posm(:,1)))';       % note, it should be transposed
cfg.fsample = 250;            % Hz
time = (1:250)/250;           % manually create a time axis
signal = sin(10*time*2*pi);   % manually create a signal
cfg.dip.signal = {repmat(signal,[length(dip_posm(:,1)) 1])};  % three trials
raw = ft_dipolesimulation(cfg); % Fieldtrip way

%% Fieldtrip Potential value, which would be used for comparison 
V1 = cell2mat(raw.trial); 