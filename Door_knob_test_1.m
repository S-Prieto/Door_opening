clear, clc

%% Load data

sim = 0;
show_images = 0;

%
if sim == 0
    folder_d = 'E:\OneDrive - Universidad de Castilla-La Mancha\Tesis_ongoing\Apertura_puertas_Andres\Picaporte\Datos\';
    % element = 'Picaporte_puerta_fino_05\Prueba_1\Prueba_1_1\Nube_Picaporte_puerta_fino_5.txt';
    element = 'Picaporte_puerta_fino_03\Prueba_1\Prueba_1_1\Nube_Picaporte_puerta_fino_03.txt';
%     element = 'Picaporte_puerta_fino_02\Prueba_1\Prueba_1_1\Nube_Picaporte_puerta_fino_02.txt';
    
else
    folder_d = 'E:\OneDrive - Universidad de Castilla-La Mancha\Tesis_ongoing\Apertura_puertas_Andres\Picaporte\Datos\Data_base\Scan_prueba\';
    element = 'Scan_600000.pcd';
end

dir = [folder_d element];

if sim == 0
    pc_raw = load(dir);
    pc_xyz = pc_raw(:,1:3);
    pc_ = pointCloud(pc_xyz);
else
    pc_raw = pcread(dir);
    pc_xyz_o = pc_raw.Location;
    pc_xyz_u = reshape(pc_xyz_o, [], 3);
    pc_xyz_u = pc_xyz_u / 10;
    pc_ = pointCloud(pc_xyz_u);
    A = [cos(-pi/2) sin(-pi/2) 0 0; -sin(-pi/2) cos(-pi/2) 0 0; 0 0 1 0; 0 0 0 1];
    tform = affine3d(A);
    pc_ = pctransform(pc_, tform);
end

figure, pcshow(pc_.Location, 'blue')
xlabel('x'), ylabel('y'), zlabel('z')

%%

max_d = 0.008; %0.005 %0.008
refV = [1 0 0]; %[1 0 0]
maxAngD = 5;
[model, inlier, outlier] = pcfitplane(pc_, max_d, refV, maxAngD);
plane = select(pc_, inlier);
[labels, numClusters] = pcsegdist(plane, 0.1);
plane_idx = mode(labels);
plane = select(plane, find(labels == plane_idx));
pc_f_aux = select(pc_, outlier);

if show_images
    figure, pcshow(plane.Location, labels)
    figure, pcshow(plane.Location,'blue')
    figure, pcshow(pc_f_aux.Location,'blue')
end

[labels, numClusters] = pcsegdist(plane, 0.002);
door_leaf = select(plane, find(labels == mode(labels)));

pc_f = pc_.Location;
idx = pc_f(:,2)<min(door_leaf.Location(:,2))...
    | pc_f(:,2)>max(door_leaf.Location(:,2));
pc_f(idx, :) = [];
pc_f = pointCloud(pc_f);

if show_images
    figure, pcshow(plane.Location, labels);
    figure, pcshow(door_leaf.Location,'blue')
    figure, pcshow(pc_f)
end

%%

interval = 0.01; %0.05 %0.01
final = 0.1; %0.05
initial = plane.Location(1,:);
point = initial;
aux = 1;
knob = [];

for i = interval*3:interval:final
    point = initial - [i 0 0];
    plane_aux_ = [[model.Normal] -dot(model.Normal, point)];
    
    model_aux = planeModel(plane_aux_);
    [pc_f, validPtCloudIndices] = removeInvalidPoints(pc_f);
    distances = evalPlane(model_aux.Parameters, pc_f.Location);
    inlier_aux = validPtCloudIndices(distances < max_d);
    plane_aux = select(pc_f, inlier_aux);
    if isempty(plane_aux.Location)
        break
    end

    [labels, numClusters] = pcsegdist(plane_aux, 0.1);
    if show_images
        figure, pcshow(plane_aux.Location, labels)
    end
    %     plane_idx = mode(labels);
    %     plane = select(plane, find(labels == 4));
    
    for j = 1:numClusters
        
       cluster_aux = select(plane_aux, find(labels == j));
       if (abs(cluster_aux.YLimits(1)-cluster_aux.YLimits(2))>0.05)
          
           knob{aux} = cluster_aux;
           
       end
        
    end
    
    aux = aux+1;
end

if length(knob)>1
    
    width = zeros(length(knob), 1);
    for i = 1:length(knob)
        
        if ~isempty(knob{i})
            width(i) = abs(knob{i}.YLimits(1)-knob{i}.YLimits(2));
        end
        
    end
    [~, knob_idx] = max(width);
    knob_slice = knob{knob_idx};
    
else
    
    knob_slice = knob{1};
    
end

if show_images
    figure, pcshow(knob_slice,'MarkerSize', 150)
end

%% KNOB PC

knob_pc = pc_f_aux.Location;
knob_pc(knob_pc(:,2)<min(knob_slice.Location(:,2)-0.02) |knob_pc(:,2) > max(knob_slice.Location(:,2))+0.02,:) = [];
knob_pc(knob_pc(:,3)<min(knob_slice.Location(:,3)-0.02) |knob_pc(:,3) > max(knob_slice.Location(:,3))+0.02,:) = [];
knob_pc = pointCloud(knob_pc);
if sim == 0
    knob_pc = pcdenoise(knob_pc);
end
[labels, numClusters] = pcsegdist(knob_pc, 0.01);
if show_images
    figure, pcshow(knob_pc.Location, 'blue');
end

%%

ptCloud = knob_pc;
wall_z_ = double(ptCloud.Location);
init_values = [min(wall_z_(:,1)) min(wall_z_(:,2)) min(wall_z_(:,3))]; 
wall_z = [wall_z_(:,1)-min(wall_z_(:,1)) wall_z_(:,2)-min(wall_z_(:,2))...
    wall_z_(:,3)-min(wall_z_(:,3))];

alto = max(wall_z(:,3));
ancho = max(wall_z(:,2));

coef = 1000;
wall_mm = [ceil(wall_z(:,1)*coef)+1 ceil(wall_z(:,2)*1000)+1 ceil(wall_z(:,3)*coef)+1];

depth_mat_front = zeros(ceil(alto*coef),ceil(ancho*coef));
% depth_norm = zeros(ceil(alto*coef),ceil(ancho*coef));
for i = 1:length(wall_z)
    depth_mat_front(wall_mm(i,3),wall_mm(i,2)) = wall_mm(i,1); 
%     depth_norm(wall_mm(i,3),wall_mm(i,1)) = wall_z_(i,5); %5
end

depth_map_front = depth_mat_front;

depth_map_front = flip(depth_map_front,1);
depth_map_front = flip(depth_map_front,2);
image_front = depth_map_front ~= 0;
% figure, imshow(image);
image_front = bwmorph(image_front, 'dilate',1);
image_front = bwmorph(image_front, 'close',1);
image_front = imfill(image_front,'holes');
% image = bwmorph(image,'erode',1);
figure, imshow(image_front);



depth_mat_top = zeros(ceil(alto*coef),ceil(ancho*coef));
% depth_norm = zeros(ceil(alto*coef),ceil(ancho*coef));
for i = 1:length(wall_z)
    depth_mat_top(wall_mm(i,1),wall_mm(i,2)) = wall_mm(i,3); 
%     depth_norm(wall_mm(i,3),wall_mm(i,1)) = wall_z_(i,5); %5
end

depth_map_top = depth_mat_top;

depth_map_top = flip(depth_map_top,1);
depth_map_top = flip(depth_map_top,2);
image_top = depth_map_top ~= 0;
% figure, imshow(image);
image_top = bwmorph(image_top, 'dilate',1);
image_top = bwmorph(image_top, 'close',1);
image_top = imfill(image_top,'holes');
% image = bwmorph(image,'erode',1);
figure, imshow(image_top);

%% ANALYTICS

p_i = max(door_leaf.Location(:,2));
p_f = min(door_leaf.Location(:,2));
door_width = abs(p_f-p_i);
k_i = max(knob_pc.Location(:,2));
k_f = min(knob_pc.Location(:,2));

if k_i < p_i  %Door handle is on the right side -> door axis on the left
%     d_h
end

door.p_i = p_i;
door.p_f = p_f;
door.width = door_width;
% door

image_top_array = sum(image_top,1);
[max_val, max_idx] = max(image_top_array);
max_idx = find(image_top_array > max_val*0.9);
axis_y_px = round((min(max_idx) + max(max_idx))/2);
z_width = find(image_front(:,axis_y_px) == 1);
axis_z_px = round((min(z_width) + max(z_width))/2);

axis_y_px_aux = size(image_top_array,2) - axis_y_px;
x_width = find(image_top(:, axis_y_px_aux) == 1);
grip_x_px = round((min(x_width) + max(x_width))/2);
if axis_y_px < size(image_top_array,2)/2
    % Door handle is on the left
    grip_y_px = round(size(image_front,2) *0.9);
    
else
    % Door handle is on the right
    grip_y_px = round(size(image_front,2) *0.1);
    
end
grip_z_px = find(image_front(:,grip_y_px) == 1, 1, 'first');

knob_y_aux = axis_y_px/size(image_front,2);
knob_width = abs(knob_pc.YLimits(1) - knob_pc.YLimits(2));
axis_y_m = max(knob_pc.YLimits) - knob_width*knob_y_aux;
knob_z_aux = axis_z_px/size(image_front,1);
knob_height = abs(knob_pc.ZLimits(1) - knob_pc.ZLimits(2));
axis_z_m = max(knob_pc.ZLimits) - knob_height*knob_z_aux;

grip_x_aux = grip_x_px/size(image_top,1);
knob_length = abs(knob_pc.XLimits(1) - knob_pc.XLimits(2));
grip_x_m = max(knob_pc.XLimits) - knob_length*grip_x_aux;
grip_y_aux = grip_y_px/size(image_front,2);
grip_y_m = max(knob_pc.YLimits) - knob_width*grip_y_aux;
grip_z_aux = grip_z_px/size(image_front,1);
grip_z_m = max(knob_pc.ZLimits) - knob_height*grip_z_aux;

figure, imshow(image_top), hold on
line([axis_y_px, axis_y_px], [1 size(image_top, 1)],'LineWidth', 2, 'Color','red')
line([1 size(image_top,2)],[grip_x_px grip_x_px],'LineWidth', 2, 'Color', 'red');
figure, imshow(image_front), hold on
line([axis_y_px, axis_y_px], [1 size(image_top, 1)],'LineWidth', 2, 'Color','red')
line([1 size(image_front,2)],[axis_z_px axis_z_px],'LineWidth', 2, 'Color', 'red')
viscircles([grip_y_px, grip_z_px], 2, 'Color', 'blue')
viscircles([axis_y_px, axis_z_px], (max(z_width) - min(z_width))/2, 'Color', 'green')

figure, pcshow(knob_pc), xlabel('x'), ylabel('y'), zlabel('z'), hold on
plot3([min(knob_pc.XLimits)*0.9 max(knob_pc.XLimits)*1.1],...
    [axis_y_m axis_y_m],[axis_z_m axis_z_m],'LineWidth', 1.5, 'Color', 'red')
scatter3(grip_x_m, grip_y_m, grip_z_m, 100,'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black')


%% REPRESENTATION

figure, pcshow(pc_f_aux.Location, 'black'), hold on
pcshow(door_leaf.Location,'blue')
pcshow(knob_pc.Location,'green','MarkerSize', 15)
plot3([min(knob_pc.XLimits)*0.9 max(knob_pc.XLimits)*1.1],...
    [axis_y_m axis_y_m],[axis_z_m axis_z_m],'LineWidth', 1.5, 'Color', 'red')
scatter3(grip_x_m, grip_y_m, grip_z_m, 100,'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black')

%%
function dis = evalPlane(model, points, varargin)
dis = abs(points * model(1:3)' + model(4));
end















