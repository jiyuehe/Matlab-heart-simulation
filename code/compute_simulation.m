function sim_u_voxel = compute_simulation(simulation_input)

create_rotor_flag = simulation_input.create_rotor_flag;
reset_time_start = simulation_input.reset_time_start;
reset_time_end = simulation_input.reset_time_end;
voxel_flag = simulation_input.voxel_flag;

% resonstruct J_stim signal
stimulus = simulation_input.stimulus;
[pacing_voxel_id,signal] = j_stim_decompress(stimulus);

% neighbor indices
geometry = simulation_input.geometry;
indices = geometry.volume.voxel_based_voxels;
indices2 = indices; % Matlab cannot index 0, so 0 is replaced with 1, 
indices2(indices2==0) = 1; % and this does not matter because the associated sign(indices(id)) will be 0

% equation parts
n_voxel = simulation_input.n_voxel;
D0 = simulation_input.D0;
d11 = zeros(n_voxel,1);
d12 = zeros(n_voxel,1);
d13 = zeros(n_voxel,1);
d21 = zeros(n_voxel,1);
d22 = zeros(n_voxel,1);
d23 = zeros(n_voxel,1);
d31 = zeros(n_voxel,1);
d32 = zeros(n_voxel,1);
d33 = zeros(n_voxel,1);
for n = 1:n_voxel
    d11(n) = D0{n}(1,1);
    d12(n) = D0{n}(1,2);
    d13(n) = D0{n}(1,3);
    d21(n) = D0{n}(2,1);
    d22(n) = D0{n}(2,2);
    d23(n) = D0{n}(2,3);
    d31(n) = D0{n}(3,1);
    d32(n) = D0{n}(3,2);
    d33(n) = D0{n}(3,3);
end
parts = zeros(n_voxel,21);
parts(:,1) = 4*d11.*sign(indices(:,1));
parts(:,2) = 4*d11.*sign(indices(:,2));
parts(:,3) = 4*d22.*sign(indices(:,3));
parts(:,4) = 4*d22.*sign(indices(:,4));
parts(:,5) = 4*d33.*sign(indices(:,5));
parts(:,6) = 4*d33.*sign(indices(:,6));
parts(:,7) = ((d11(indices2(:,1))-d11(indices2(:,2))).*sign(indices(:,1)).*sign(indices(:,2)) + (d21(indices2(:,3))-d21(indices2(:,4))).*sign(indices(:,3)).*sign(indices(:,4)) + (d31(indices2(:,5))-d31(indices2(:,6))).*sign(indices(:,5)).*sign(indices(:,6))) .* sign(indices(:,1)).*sign(indices(:,2));
parts(:,8) = ((d12(indices2(:,1))-d12(indices2(:,2))).*sign(indices(:,1)).*sign(indices(:,2)) + (d22(indices2(:,3))-d22(indices2(:,4))).*sign(indices(:,3)).*sign(indices(:,4)) + (d32(indices2(:,5))-d32(indices2(:,6))).*sign(indices(:,5)).*sign(indices(:,6))) .* sign(indices(:,3)).*sign(indices(:,4));
parts(:,9) = ((d13(indices2(:,1))-d13(indices2(:,2))).*sign(indices(:,1)).*sign(indices(:,2)) + (d23(indices2(:,3))-d23(indices2(:,4))).*sign(indices(:,3)).*sign(indices(:,4)) + (d33(indices2(:,5))-d33(indices2(:,6))).*sign(indices(:,5)).*sign(indices(:,6))) .* sign(indices(:,5)).*sign(indices(:,6));
parts(:,10) = 2*d12.*sign(indices(:,7)).*sign(indices(:,9));
parts(:,11) = 2*d12.*sign(indices(:,8)).*sign(indices(:,10));
parts(:,12) = 2*d13.*sign(indices(:,15)).*sign(indices(:,17));
parts(:,13) = 2*d13.*sign(indices(:,16)).*sign(indices(:,18));
parts(:,14) = 2*d23.*sign(indices(:,11)).*sign(indices(:,13));
parts(:,15) = 2*d23.*sign(indices(:,12)).*sign(indices(:,14));
parts(:,16) = simulation_input.tau_open_voxel;
parts(:,17) = simulation_input.tau_close_voxel;
parts(:,18) = simulation_input.tau_in_voxel;
parts(:,19) = simulation_input.tau_out_voxel;
parts(:,20) = simulation_input.v_gate_voxel;
parts(:,21) = simulation_input.c_voxel;

% compute simulation
u_current = zeros(n_voxel,1); % initial value 0, set all voxel at rest
h_current = ones(n_voxel,1); % initial value 1, set all voxel at rest
h_next = zeros(n_voxel,1);
t_final = simulation_input.t_final; % unit: ms
dt = simulation_input.dt;
T = t_final/dt; % number of simulation time steps
delta = geometry.volume.delta; % voxel spatial distance
id_save = 1;
sim_u_voxel = zeros(n_voxel,t_final); % sampling frequency at 1 kHz
for t = 1:T-1
    do_flag = 1;    
    if do_flag == 1 && mod(t,T/5) == 1
        disp(['simulation ',num2str((t-1)/T*100),'%']);
    end    

    % compute diffusion term
    diffusion_term = 1/(4*delta^2) * ...
        (parts(:,1) .* (u_current(indices2(:,1))-u_current) + parts(:,2) .* (u_current(indices2(:,2))-u_current) + ...
        parts(:,3) .* (u_current(indices2(:,3))-u_current) + parts(:,4) .* (u_current(indices2(:,4))-u_current) + ...
        parts(:,5) .* (u_current(indices2(:,5))-u_current) + parts(:,6) .* (u_current(indices2(:,6))-u_current) + ...
        parts(:,7) .* (u_current(indices2(:,1))-u_current(indices2(:,2))) + ...
        parts(:,8) .* (u_current(indices2(:,3))-u_current(indices2(:,4))) + ...
        parts(:,9) .* (u_current(indices2(:,5))-u_current(indices2(:,6))) + ...
        parts(:,10) .* (u_current(indices2(:,7))-u_current(indices2(:,9))) + parts(:,11) .* (u_current(indices2(:,10))-u_current(indices2(:,8))) + ...
        parts(:,12) .* (u_current(indices2(:,15))-u_current(indices2(:,17))) + parts(:,13) .* (u_current(indices2(:,18))-u_current(indices2(:,16))) + ...
        parts(:,14) .* (u_current(indices2(:,11))-u_current(indices2(:,13))) + parts(:,15) .* (u_current(indices2(:,14))-u_current(indices2(:,12))) );
    diffusion_term = parts(:,21) .* diffusion_term;
    
    if create_rotor_flag == 1 && t >= reset_time_start && t <= reset_time_end
        % set region 4 to at rest
        vol_id = (voxel_flag==4);
        u_current(vol_id) = 0;
        h_current(vol_id) = 1;
    end
    
    % compute the next time step value of u
    J_stim = zeros(n_voxel,1);
    J_stim(pacing_voxel_id) = signal(:,t);
    u_next = ((h_current .* u_current.^2 .* (1-u_current) ./ parts(:,18)) + ...
        (-u_current ./ parts(:,19)) + J_stim + diffusion_term) * dt + u_current;
    
    % compute the next time step value of h
    h_next_1 = ((1-h_current) ./ parts(:,16)) * dt + h_current;
    h_next_2 = (-h_current ./ parts(:,17)) * dt + h_current;
    id_1 = u_current < parts(:,20);
    id_2 = u_current >= parts(:,20);
    h_next(id_1) = h_next_1(id_1);
    h_next(id_2) = h_next_2(id_2);

    % update value
    u_current = u_next;
    h_current = h_next;

    % save value
    if mod(t,1/dt) == 1
        sim_u_voxel(:,id_save) = u_current;
        id_save = id_save + 1;
    end
end

end
