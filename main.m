clear all;
close all;
clc;

%% Important References
% http://www.fieldtriptoolbox.org/example/compute_forward_simulated_data_and_apply_a_dipole_fit/
% http://www.fieldtriptoolbox.org/reference/ft_prepare_leadfield/
% http://www.fieldtriptoolbox.org/reference/ft_prepare_sourcemodel/

%% Add path of fieldtrip toolbox
addpath('/Users/aarushigupta/Desktop/BEM_Amita/Real Head Model/fieldtrip-20190224/fieldtrip-20190224')
ft_defaults  % initiate fieldtrip

%% STEP-1 Reading MRI file - anatomical data
mri = ft_read_mri('standard_mri.mat');  % E:\EEG_DATA\fieldtrip-20190224\fieldtrip-20190224\template\headmodel\standard_mri.mat 
cfg        = [];    % To plot MRI define configuration
cfg.method = 'ortho'; % cfg.method could be 'slice', 'ortho', 'glassbrain'
cfg.colorbar = 'no';
ft_sourceplot(cfg, mri);  
ft_determine_coordsys(mri, 'interactive', 'no');

%% STEP-2 Segmentation of the anatomical information into three different tissue types
cfg        = [];    % cfg: configuration of the function that was used to create vol
cfg.output       = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);
% save segmentedmri segmentedmri

%% STEP-3 Define number of vertices of each segmented tissue to prepare triangulation mesh 
load segmentedmri
cfg        = [];
cfg.tissue = {'brain','skull','scalp'};
cfg.numvertices = [1500 1500 1500]; %1500 vertices over each surface
bnd = ft_prepare_mesh(cfg,segmentedmri); % The bnd field contains information about the mesh
% save bnd bnd

%% Prepare Headmodel - bemcp method i.e Boundary element method
load bnd
cfg = [];
cfg.method = 'bemcp'; % You can also specify 'openmeeg', 'dipoli', or another method.
cfg.conductivity = [0.3300 0.0041 0.3300];
vol = ft_prepare_headmodel(cfg,bnd);
% save vol vol

%% Visualization plot mesh
load vol
figure
ft_plot_mesh(vol.bnd(1), 'facecolor',[0 0 1], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'facecolor',[0 1 1],'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facealpha', 0.3,'facecolor',[0.4 0.6 0.4]);

%% Electrode placement over the scalp surface
elec = ft_read_sens('standard_1020.elc');  % for electrodes - 97 electrodes - default - 'standard_1020.elc'
figure;
ft_plot_mesh(vol.bnd(3), 'edgecolor','none','facealpha',0.3,'facecolor',[0.4 0.6 0.4]);
hold on;
ft_plot_sens(elec);     % optional properties - style,marker

