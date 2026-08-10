nodes = load('nodes.dat');        % [x y z id]
elements = load('tetraedre.dat');  % [n1 n2 n3 n4]

% coords and indice
coords = nodes(:,1:3);
node_ids = nodes(:,4);

% the sphere
center = [0, 0, 0];
radius = 0.5;

% nodes in sphere
distances = sqrt(sum((coords - center).^2, 2));
inside = distances < radius;

% new nodes 
coords_new = coords(~inside, :);
node_ids_new = node_ids(~inside);

% new tetra
keep = ~any(ismember(elements, node_ids(inside)), 2);
elements_new = elements(keep, :);

% new number
[~, loc] = ismember(elements_new, node_ids_new);
elements_reindexed = loc;


fprintf('number of nodes: %d\n', size(coords,1));
fprintf('number of new nodes: %d\n', size(coords_new,1));
fprintf('number of tetra: %d\n', size(elements,1));
fprintf('number of new tetra: %d\n', size(elements_reindexed,1));

% plot
figure;
tetramesh(elements_reindexed, coords_new);
axis equal;
title('Mesh after removing inner sphere');


save('remaining_nodes.dat','coords_new','-ascii');
save('remaining_elements.dat','elements_reindexed','-ascii');