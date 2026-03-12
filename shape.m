% Load data from file
fid = fopen('original/nodes10.dat');
data = fscanf(fid, '%g %g %g', [3, Inf]);
A = data';
fclose(fid);

% Define radius
R = 0.5;

% Extract points outside the sphere of radius R
M = [];

for i = 1:length(A)

    if norm(A(i, :)) > R
        M = [M; A(i, :)];
    end

end

% Define P = points in M plus a sphere of points at radius R
[x1, y1, z1] = sphere(12);
Psphere = [R * x1(:) R * y1(:) R * z1(:)];
Psphere = unique(Psphere, 'rows');
P = [M; Psphere];

fprintf('%d points in tetra, %d points in tetra with hole, %d points in tetra with hole plus spherical cavity.\n', ...
    length(A), length(M), length(P));

% Create alpha shape
shp = alphaShape(P(:, 1), P(:, 2), P(:, 3), R * 0.8);
plot(shp)

% Get triangulation
[tri, loc] = alphaTriangulation(shp);

% Optional: save triangulation and nodes
% fid = fopen('nodess.dat','w');
% fprintf(fid,'%i\t %i\t %i\n',loc);
% fclose(fid);
% fid = fopen('tetreades.dat','w');
% fprintf(fid,'%i\t %i\t %i\t %i\n',tri);
% fclose(fid);

% Report mesh statistics
numtetrahedra = size(tri)
numnodes = size(loc)

% Plot shape with axis labels
%plot(shp)
%axis equal
%xlabel('X'); ylabel('Y'); zlabel('Z')
%view(90, -10) % rotate view of shape

% Create PDE model and mesh
model = createpde;
[G, mesh] = geometryFromMesh(model, loc', tri');

% Save results
save('code/tetra.dat', 'tri', '-ascii')
save('code/coordinates.dat', 'loc', '-ascii')

axis equal


% ── Update header file with mesh statistics ───────────────────────────────────
gm = numtetrahedra(1);
nt = numnodes(1);
update_header(gm, nt);

% ─────────────────────────────────────────────────────────────────────────────
function update_header(gm, nt)
    HEADER_PATH = 'code/fondements_keltoum.h';

    fid = fopen(HEADER_PATH, 'r');
    content = fread(fid, '*char')';
    fclose(fid);

    replacement = sprintf('gm= %d, nt= %d', gm, nt);
    pattern = '(gm=\s*)\d+(,\s*nt=\s*)\d+';
    new_content = regexprep(content, pattern, replacement, 'once');

    fid = fopen(HEADER_PATH, 'w');
    fwrite(fid, new_content);
    fclose(fid);

    fprintf('Updated %s: gm=%d, nt=%d\n', HEADER_PATH, gm, nt);
end
