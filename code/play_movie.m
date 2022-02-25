function play_movie(directory,file_name,movie_type,data,data_min,data_max,frame_id)

home_dir = directory.home_dir;
result_dir = directory.result_dir;

signal = data.signal;
az = data.az; % view angle
el = data.el; % view angle

if length(frame_id) > 1 
    cd(result_dir);
    writerObj = VideoWriter(file_name);
    writerObj.FrameRate = 10; % frames per second
    open(writerObj);
end
figure('units','pixels','position',[100 100 1280 720]);
movegui('center');
set(gcf,'color','w');
if length(frame_id) > 1 
    pause(2); % pause a little bit to wait for the figure to fully expanded
end
set(gcf,'Renderer','zbuffer');

if strcmp(movie_type,'point_cloud')
    vertex = data.voxel;
    h = scatter3(vertex(:,1),vertex(:,2),vertex(:,3),10,'b','filled'); 
elseif strcmp(movie_type,'triangular_mesh')
    vertex = data.atrium.vertex;
    face = data.atrium.face;
    h = patch('Faces',face,'Vertices',vertex,'FaceColor','interp','EdgeColor',0.3*[1 1 1]);
elseif strcmp(movie_type,'triangular_mesh_phase_singularity')
    vertex = data.atrium.vertex;
    face = data.atrium.face;
    faceType = data.faceType;    
    h1 = patch('Vertices',vertex,'Faces',face,'FaceColor','interp','FaceColor','w','EdgeColor','k','edgealpha',0.2);
    h2 = patch('Vertices',vertex,'Faces',face,'FaceColor','flat','FaceColor','w','EdgeColor','k','edgealpha',0.2);
end
% hold on;
% pacing_vertex_id = find(data.vertex_flag==2);
% scatter3(vertex(pacing_vertex_id,1),vertex(pacing_vertex_id,2),vertex(pacing_vertex_id,3),...
%     30,'r','filled'); % mark pacing vertices
% hold off;
view(az,el);
axis tight equal vis3d off;
rotate3d on;

for n = frame_id
    [~,movie_color] = convert_data_to_color(data_max,data_min,data_min,signal(n,:));
    
    if strcmp(movie_type,'point_cloud')
        set(h,'CData',movie_color);
    elseif strcmp(movie_type,'triangular_mesh')
        % set back face black
%         back_vertex_id = find_vertices_id_on_the_back(face,vertex);
%         movie_color(back_vertex_id,:) = 0;        
        
        set(h,'FaceVertexCData',movie_color);
    elseif strcmp(movie_type,'triangular_mesh_phase_singularity')
        set(h1,'FaceColor','interp','FaceVertexCData',movie_color);
        
        faceID = [];
        facesColor = [];
        
        lineFaces = faceType{n}.rg; % wavefront line
        faceID(end+1:end+length(lineFaces)) = lineFaces;
        facesColor(end+1:end+length(lineFaces),:) = repmat([1 1 0],length(lineFaces),1); % yellow
        
        cwFaces = faceType{n}.cw; % clockwise phase singularity
        faceID(end+1:end+length(cwFaces)) = cwFaces;
        facesColor(end+1:end+length(cwFaces),:) = repmat([1 0 0],length(cwFaces),1); % red
        
        ccwFaces = faceType{n}.ccw; % counterclockwise phase singularity
        faceID(end+1:end+length(ccwFaces)) = ccwFaces;
        facesColor(end+1:end+length(ccwFaces),:) = repmat([0 0 1],length(ccwFaces),1); % blue
        
        set(h2,'FaceColor','flat','Faces',face(faceID,:),'FaceVertexCData',facesColor);
    end    
    
    title(['time: ',num2str(n-1),' ms']);
    
    if length(frame_id) > 1 
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
end

if length(frame_id) > 1 
    close(writerObj);
    close;
    cd(home_dir);
end
    
end
    