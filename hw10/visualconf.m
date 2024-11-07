function visualconf(data,fig)
Na = size(data,2);
rstar = 0.1;
pad = 0.05 + rstar;
cr = 0.02; % radius of cylinder
dtol = 1e-2;
r = 2^(1/6);
r2 = r*sqrt(2);
tol = 0.1*r;

t = linspace(-1,1,50);
[X, Y, Z] = meshgrid(t,t,t);
V = sqrt(X.^2 + Y.^2 + Z.^2);
[~,verts] = isosurface(X,Y,Z,V,rstar);
nverts = size(verts,1);
e = ones(nverts,1);

% colors
vcol = [0 0 0]/255; % color for vertices
ecol = [32 32 32]/255; % color for edges
fcol = [255 0 0]/255; % color for faces


t = linspace(-0.5,0.5,50);
[X, Y, Z] = meshgrid(t,t,t);
V = sqrt(X.^2 + Y.^2 + Z.^2);
[bfaces,bverts] = isosurface(X,Y,Z,V,rstar); % faces and vertices of balls
nbverts = size(bverts,1);
eb = ones(nbverts,1);
V = sqrt(X.^2 + Y.^2);
[cfaces,cverts] = isosurface(X,Y,Z,V,cr); % faces and vertices of balls
ncverts = size(cverts,1);
ec = ones(ncverts,1); 

%% find polyhedron
x = data;
cx = sum(x,2)/Na;
x = x - cx*ones(1,Na);

%% find the distance matrix
for j = 1 : Na
    dist(j,:) = sqrt(sum((x - x(:,j)*ones(1,Na)).^2,1));
end

%% find faces
tfaces = 0;
sfaces = 0;
for j=1: Na-2
    dd = dist(j,j+1:Na);
    ind1 = find(abs(dd - r) < tol);
    if ~isempty(ind1)
        li = length(ind1);
        % find equilateral triangle faces
        for k = 1:li
            i1 = ind1(k);
            dd1 = dist(j,j+i1+1:Na);
            dd2 = dist(j+i1,j+i1+1:Na);
            % find equilateral triangle faces
            ind = find(abs(dd1 - r) < tol & abs(dd2 - r) < tol);
            if ~isempty(ind)
                l=length(ind);
                for i=1:l
                    tfaces = tfaces+1;
                    tfac(tfaces,:) = [j,j+i1,j+i1+ind(i)];
                end
            end
        end
        % find square faces
        if li >= 2
            for k1 = 1:li
                i1 = ind1(k1);
                dd1 = dist(j+i1,j+1:Na);
                for k2 = 1:li
                    i2 = ind1(k2);
                    if abs(dist(j+i1,j+i2) - r2) < tol
                        dd2 = dist(j+i2,j+1:Na);
                        ind = find(abs(dd1 - r) < tol & abs(dd2 - r) < tol ...
                            & abs(dd - r2) < tol );
                        if ~isempty(ind)
                            l = length(ind);
                            for i=1:l
                                sfaces = sfaces+1;
                                sfac(sfaces,:) = [j,j+i1,j+ind(i),j+i2];
                            end
                        end
                    end
                end
            end
        end
    end
end



v = x';
%% visualize
figure(fig); clf;
hold on;

ht = patch('vertices',v,'faces',tfac,'facecolor',fcol,'EdgeColor','none','LineStyle','none');
alpha(ht,0.35); %0.35
hold on
% if sfaces > 0
%     hs = patch('vertices',v,'faces',sfac,'facecolor',fcol,'EdgeColor','none','LineStyle','none');
%     alpha(hs,0.35) %0.35
% end


%% draw vertices and edges
for j = 1 : Na 
    xyz = x(:,j);
    hv = patch('Vertices',bverts + eb*xyz','Faces',bfaces,'Facecolor',vcol,'EdgeColor','none','LineStyle','none');
end
for i = 1 : Na - 1
    for j = i + 1 : Na
        if abs(norm(x(:,i) - x(:,j)) - r) < tol
            x0 = x(:,i);
            x1 = x(:,j);
            u = x1 - x0;
            u = u/norm(u); % unit vector parallel to the the vector x0 --> x1
            R = make_rotation_matrix(u); % makes rotation matrix 
            % for turn around axis normal to Oz cross u by angle between (Oz,u)
            cv = cverts*R';
            cv = cv + 0.5*ec*(x0+x1)'; % shift start of the cylinder to 0.5*(x0 + x1);
            patch('Vertices',cv,'Faces',cfaces,'Facecolor',ecol,'EdgeColor','none','LineStyle','none');
        end
    end
end
view(3)
set(gca,'DataAspectRatio',[1,1,1]);
camlight('headlights')
lighting gouraud
axis off
drawnow;
end
%%
%%
function R = make_rotation_matrix(u)
I = eye(3);
e = [0;0;1];
co = u'*e;
u = cross(e,u);
si = norm(u);
if si < 1e-12
    if co > 0
        R = I;
    else
        R = I;
        R(3,3) = -1;
    end
    return
end
u = u/si;    
ucr = [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
uu = u*u';
R = I*co + ucr*si + (1 - co)*uu;
end
