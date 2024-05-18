clear; close all; clc;
addpath(genpath('code'));

directory.home_dir = pwd;
directory.data_dir = fullfile(directory.home_dir,'data');
directory.result_dir = fullfile(directory.home_dir,'result');

% setting
geometry_id = 1; % 1: atrium, 2: slab

% load data
cd(directory.data_dir);
if geometry_id == 1
    load('atrium_simulation_input.mat');
    load('atrium_geometry.mat');
elseif geometry_id == 2
    load('slab_simulation_input.mat');
    load('slab_geometry.mat');
end
cd(directory.home_dir);
t_final = 500; % unit: ms
dt = 0.01; % unit: ms

% run simulation
simulation_result_file_name = 'simulation_result';
tic; % start timer
sim_v_voxel = compute_simulation(simulation_input);
action_potential_voxel = sim_v_voxel';
toc; % end timer, will show computation time in command window

% show simulation result
n_voxel = size(action_potential_voxel,2);
v_gate_voxel = ones(n_voxel,1) * 0.13;

action_potential_voxel_phase = zeros(size(action_potential_voxel));
for id = 1:n_voxel
    action_potential_voxel_phase(:,id) = create_movie_phase(action_potential_voxel(:,id),v_gate_voxel(id));
end

% movie played on point cloud
debug_plot = 1;
if debug_plot == 1
    data.voxel = geometry.volume.voxel;
    data.signal = action_potential_voxel_phase;
    data.az = 10;
    data.el = 60;
    play_movie(directory,['simulation_point_cloud_',num2str(geometry_id)],'point_cloud',data,v_gate_voxel(1),1,1:10:size(action_potential_voxel_phase,1));
end

% movie played on triangular mesh
debug_plot = 1;
if debug_plot == 1
    sim_v_mesh = sim_v_voxel(geometry.id_mapping.voxel_id_for_each_vertex,:);   
    action_potential = sim_v_mesh';
    
    action_potential_phase = zeros(size(action_potential));
    v_gate_vertex = v_gate_voxel(geometry.id_mapping.voxel_id_for_each_vertex);
    for id = 1:size(action_potential,2)
        action_potential_phase(:,id) = create_movie_phase(action_potential(:,id),v_gate_vertex(id));
    end
    
    data.atrium = geometry.atrium;
    data.signal = action_potential_phase;
    data.az = 10;
    data.el = 60;
    play_movie(directory,['simulation_triangular_mesh_',num2str(geometry_id)],'triangular_mesh',data,v_gate_vertex(1),1,1:10:size(action_potential_phase,1));
end
