% Load data from file
fid = fopen('original/nodes10.dat', 'r');
data = fscanf(fid, '%f %f %f', [3, Inf]);
A = data';
fclose(fid);

% Sphere center and radius
Center = [0 0 0];
R = 1.5;

% Remove points strictly inside the sphere
distances = sqrt(sum((A - Center) .^ 2, 2));
RemainingPoints = A(distances > R, :);

fprintf('Number of points before removal: %d\n', size(A, 1));
fprintf('Number of points after removal: %d\n', size(RemainingPoints, 1));

% Generate tetrahedral mesh (nodes + elements)
DT = delaunayTriangulation(RemainingPoints);
nodes = DT.Points;
elements = DT.ConnectivityList;

% Report mesh statistics
numtetrahedra = size(elements);
numnodes = size(nodes);

% Optional plotting (comment out for batch runs)
% figure
% tetramesh(DT, 'FaceAlpha', 0.2)
% axis equal
% title('Mesh after removing interior points')
%
% shp = alphaShape(RemainingPoints, 2);
% figure
% plot(shp, 'FaceColor', 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
% axis equal
% title('3D shape with internal cavity')

% Save results in the expected locations
save('code/tetra.dat', 'elements', '-ascii')
save('code/coordinates.dat', 'nodes', '-ascii')


% ── Update header file with mesh statistics ───────────────────────────────────
nt = numtetrahedra(1);
gm = numnodes(1);
update_header(gm, nt);

% ─────────────────────────────────────────────────────────────────────────────
function update_header(gm, nt)
    fid = fopen('code/fondements_keltoum.h');
    content = fread(fid, '*char')';
    fclose(fid);

    replacement = sprintf('gm= %d, nt= %d', gm, nt);
    pattern = '(gm=\s*)\d+(,\s*nt=\s*)\d+';
    new_content = regexprep(content, pattern, replacement, 'once');

    fid = fopen('code/fondements_keltoum.h', 'w');
    fwrite(fid, new_content);
    fclose(fid);

    fprintf('Updated %s: gm=%d, nt=%d\n', 'fondements_keltoum.h', gm, nt);
end
