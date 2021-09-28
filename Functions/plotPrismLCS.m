function [] = plotPrismLCS(pc,LCS,tcoef,index,h_star,sigma_pcd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

count = pc.Count;

if(not(isempty(pc.Color)))
    colorlcs(1,:) = pc.Color(index,:);
    for j=1:count
        if j<index
            colorlcs(j+1,:) = pc.Color(j,:);
        end
        if j>index
            colorlcs(j,:) = pc.Color(j,:);
        end
    end
    pcshow(LCS,colorlcs);
else
    pcshow(LCS);
end

for q=1:4
    
    switch q
        case 1
            h = h_star(1);

            vertices_matrix = [0 0 -max(h/2,3*sigma_pcd); h 0 -max(h/2,3*sigma_pcd); h h -max(h/2,3*sigma_pcd); 0 h -max(h/2,3*sigma_pcd); ...
                                            0 0 max(h/2,3*sigma_pcd); h 0 max(h/2,3*sigma_pcd); h h max(h/2,3*sigma_pcd); 0 h max(h/2,3*sigma_pcd)];

            prism = zeros(count,3,8);

            for i=1:8
                coord = vertices_matrix(i,:);
                temp = coord/tcoef;
                point = pc.Location(index,:);
                prism(index,:,i) = temp;
            end
            
            vertices = [prism(index,:,1); prism(index,:,2); prism(index,:,3); prism(index,:,4); prism(index,:,5); prism(index,:,6); prism(index,:,7); prism(index,:,8)]; 
            faces_matrix = [1 2 6 5; 2 3 7 6; 3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
            h = patch('Vertices',vertices,'Faces',faces_matrix,...
            'FaceVertexCData',hsv(8), ...
            'FaceAlpha',0.5, ...
            'FaceColor','flat');
            view(3)
            xlabel('x', 'FontName', 'courier new', 'FontWeight', 'bold');
            ylabel('y', 'FontName', 'courier new', 'FontWeight', 'bold');
            zlabel('z', 'FontName', 'courier new', 'FontWeight', 'bold');
            hold on;
            grid on;
            box on;
    
        case 2
            h = h_star(2);

            vertices_matrix = [0 0 -max(h/2,3*sigma_pcd); -h 0 -max(h/2,3*sigma_pcd); -h h -max(h/2,3*sigma_pcd); 0 h -max(h/2,3*sigma_pcd); ...
                                            0 0 max(h/2,3*sigma_pcd); -h 0 max(h/2,3*sigma_pcd); -h h max(h/2,3*sigma_pcd); 0 h max(h/2,3*sigma_pcd)];

            prism = zeros(count,3,8);

            for i=1:8
                coord = vertices_matrix(i,:);
                temp = coord/tcoef;
                point = pc.Location(index,:);
                prism(index,:,i) = temp;
            end
            
            vertices = [prism(index,:,1); prism(index,:,2); prism(index,:,3); prism(index,:,4); prism(index,:,5); prism(index,:,6); prism(index,:,7); prism(index,:,8)]; 
            faces_matrix = [1 2 6 5; 2 3 7 6; 3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
            h = patch('Vertices',vertices,'Faces',faces_matrix,...
            'FaceVertexCData',hsv(8), ...
            'FaceAlpha',0.5, ...
            'FaceColor','flat');
            view(3)
            xlabel('x', 'FontName', 'courier new', 'FontWeight', 'bold');
            ylabel('y', 'FontName', 'courier new', 'FontWeight', 'bold');
            zlabel('z', 'FontName', 'courier new', 'FontWeight', 'bold');
            hold on;
            grid on;
            box on;
            
        case 3
            h = h_star(3);

            vertices_matrix = [0 0 -max(h/2,3*sigma_pcd); -h 0 -max(h/2,3*sigma_pcd); -h -h -max(h/2,3*sigma_pcd); 0 -h -max(h/2,3*sigma_pcd); ...
                                         0 0 max(h/2,3*sigma_pcd); -h 0 max(h/2,3*sigma_pcd); -h -h max(h/2,3*sigma_pcd); 0 -h max(h/2,3*sigma_pcd)];

            prism = zeros(count,3,8);

            for i=1:8
                coord = vertices_matrix(i,:);
                temp = coord/tcoef;
                point = pc.Location(index,:);
                prism(index,:,i) = temp;
            end
            
            vertices = [prism(index,:,1); prism(index,:,2); prism(index,:,3); prism(index,:,4); prism(index,:,5); prism(index,:,6); prism(index,:,7); prism(index,:,8)]; 
            faces_matrix = [1 2 6 5; 2 3 7 6; 3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
            h = patch('Vertices',vertices,'Faces',faces_matrix,...
            'FaceVertexCData',hsv(8), ...
            'FaceAlpha',0.5, ...
            'FaceColor','flat');
            view(3)
            xlabel('x', 'FontName', 'courier new', 'FontWeight', 'bold');
            ylabel('y', 'FontName', 'courier new', 'FontWeight', 'bold');
            zlabel('z', 'FontName', 'courier new', 'FontWeight', 'bold');
            hold on;
            grid on;
            box on;
            
        case 4
            h = h_star(4);

            vertices_matrix = [0 0 -max(h/2,3*sigma_pcd); h 0 -max(h/2,3*sigma_pcd); h -h -max(h/2,3*sigma_pcd); 0 -h -max(h/2,3*sigma_pcd); ...
                                            0 0 max(h/2,3*sigma_pcd); h 0 max(h/2,3*sigma_pcd); h -h max(h/2,3*sigma_pcd); 0 -h max(h/2,3*sigma_pcd)];

            prism = zeros(count,3,8);

            for i=1:8
                coord = vertices_matrix(i,:);
                temp = coord/tcoef;
                point = pc.Location(index,:);
                prism(index,:,i) = temp;
            end
            
            vertices = [prism(index,:,1); prism(index,:,2); prism(index,:,3); prism(index,:,4); prism(index,:,5); prism(index,:,6); prism(index,:,7); prism(index,:,8)]; 
            faces_matrix = [1 2 6 5; 2 3 7 6; 3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
            h = patch('Vertices',vertices,'Faces',faces_matrix,...
            'FaceVertexCData',hsv(8), ...
            'FaceAlpha',0.5, ...
            'FaceColor','flat');
            view(3)
            xlabel('x', 'FontName', 'courier new', 'FontWeight', 'bold');
            ylabel('y', 'FontName', 'courier new', 'FontWeight', 'bold');
            zlabel('z', 'FontName', 'courier new', 'FontWeight', 'bold');
            hold on;
            grid on;
            box on;
    end
end

end

