function LJ_line_search()
fsz = 20; % fontsize
Na = 7; % the number of atoms
rstar = 2^(1/6); % argument of the minimum of the Lennard-Jones pair potential V(r) = r^(-12) - r^(-6)
tol = 1e-9; % stop iterations when || grad f|| < tol
iter_max = 500; % the maximal number of iterations
draw_flag = 1; % if draw_flag = 1, draw configuration at every iteration
% parameters for backtracking line search
c = 0.1;
rho = 0.9;
%% Set up the initial configuration

% Four lical minima of LJ7:
% f1 = -16.50538417 Pentagonal bipyramid 
% f2 = -15.93504306 Capped octahedron 
% f3 = -15.59321094 Tricapped tetrahedron 
% f4 = -15.53306005 Bicapped trigonal bipyramid

% Options: model = 0,1,2,3, or 4.
% Model 0 corresponds to a random initialization.
% Models 1--4 set the system up close to the corresponding local minima
% listed above.

model = 1;
if model > 0
    Na = 7;
end
xyz = initial_configuration(model,Na,rstar);
drawconf(xyz,1);

x = remove_rotations_translations(xyz);
LJhess(x)
%% start minimization
% choose algorithm
% direction = 1: steepest descent
% direction = 2: Newton
direction = 2;

f = LJpot(x);
g = LJgrad(x);
norm_g = norm(g);
fprintf("Initially, f = %d, ||grad f|| = %d\n",f,norm_g);
iter = 1;

fvals = zeros(iter_max);
fvals(1) = f;
ngvals = zeros(iter_max);
ngvals(1) = norm_g;

fail_flag = 0;

while norm_g > tol && iter < iter_max 
    % choose search direction
    switch direction
        case 1 % steepest descent
            p = -g;
            dir = "SD";
        case 2 % Newton
            H = LJhess(x);
            [~,flag] = chol(H);
            if flag == 0 % H is SPD, use Newton's direction
                p = -H\g; 
                dir = "Newton";
            else % use the steepest descent direction
                p = -g;
                dir = "SD";
            end
        otherwise
            return
    end
    % normalize the search direction if its length greater than 1
    norm_p = norm(p);
    if norm_p > 1
        p = p/norm_p;
    end
    % do backtracking line search along the direction p
    a = 1;
    f_temp = LJpot(x + a*p);
    cpg = c*p'*g;
    while f_temp > f + a*cpg % check Wolfe's condition 1
        a = a*rho;
        if a < 1e-14
            fprintf("line search failed\n");
            iter = iter_max;
            fail_flag = 1;
            break;
        end
        f_temp = LJpot(x + a*p);        
    end
    x = x + a*p;
    f = LJpot(x);
    g = LJgrad(x);
    norm_g = norm(g);
    fprintf("iter %d : dir = %s, f = %d, ||grad f|| = %d, step length = %d\n",iter,dir,f,norm_g,a);
    if draw_flag == 1
        % restore all coordinates
        xyz = reshape([0;0;0;x(1);0;0;x(2:3);0;x(4:end)],3,Na);
        drawconf(xyz,1);
    end
    iter = iter + 1;
    fvals(iter) = f;
    ngvals(iter) = norm_g;
    if fail_flag == 1
        break;
    end
end
xyz = reshape([0;0;0;x(1);0;0;x(2:3);0;x(4:end)],3,Na);
drawconf(xyz,1);

figure(2);
clf;
hold on;
grid on;
plot(0:(iter-1),fvals(1:iter),'Linewidth',2)
set(gca,'Fontsize',fsz);
xlabel('Iteration #','FontSize',fsz);
ylabel('Function values','FontSize',fsz);

figure(3);
clf;
hold on;
grid on;
plot(0:(iter-1),ngvals(1:iter),'Linewidth',2)
set(gca,'Yscale','log','Fontsize',fsz);
xlabel('Iteration #','FontSize',fsz);
ylabel('||grad f||','FontSize',fsz);

if fail_flag == 0
    xyz = reshape([0;0;0;x(1);0;0;x(2:3);0;x(4:end)],3,Na);
    visualconf(xyz,4);
end

figure(5);
clf;
hold on;
grid on;
evals = sort(eig(LJhess(x)),'ascend');
plot(evals,'.','Markersize',20)
if evals(1) > 0
    set(gca,'Yscale','log','Fontsize',fsz);
else
    set(gca,'Yscale','linear','Fontsize',fsz);
end    
xlabel('index','FontSize',fsz);
ylabel('Eigenvalues of the Hessian','FontSize',fsz);

end

%% make initial configuration
function xyz = initial_configuration(model,Na,rstar)
xyz = zeros(3,Na);
switch(model)
    case 1 % Pentagonal bipyramid
        p5 = 0.4*pi;
        he = sqrt(1 - (0.5/sin(0.5*p5))^2);
        for k = 1 : 5
            xyz(1,k) = cos((k-1)*p5); 
            xyz(2,k) = sin((k-1)*p5); 
            xyz(3,k) = 0;  
        end
        xyz(3,6) = he;
        xyz(3,7) = -he;

case 2 % Capped octahedron
        r = 1/sqrt(2);
        p4 = 0.5*pi;
        pp = p4/2;
        p0 = pi*1.5 - pp;
        x0 = sin(pp);
        y0 = cos(pp);
        z0 = 0;
        for k = 1 : 4
            xyz(1,k) = x0 + r*cos(p0 + (k-1)*p4);
            xyz(2,k) = y0 + r*sin(p0 + (k-1)*p4);
            xyz(3,k) = z0;
        end
        xyz(:,5) = [x0, y0, z0 + r]';
        xyz(:,6) = [x0, y0, z0 - r]';
        xyz(:,7) = [3*x0, y0, z0 + r]';

    case 3  % Tricapped tetrahedron
        p3 = 2*pi/3;
        pp = p3/2;
        r = 1/sqrt(3);
        beta = 0.5*pi -asin(1/3) - acos(1/sqrt(3));
        r1 = cos(beta);
        p0 = 1.5*pi - pp;
        x0 = sin(pp);
        y0 = cos(pp);
        z0 = 0;
        for k = 1 : 3
            xyz(1,k) = x0 + r*cos(p0 + (k-1)*p3);
            xyz(2,k) = y0 + r*sin(p0 + (k-1)*p3);
            xyz(3,k) = z0;
            xyz(1,k + 3) = x0 + r1*cos(p0 + pp + (k-1)*p3);
            xyz(2,k + 3) = y0 + r1*sin(p0 + pp + (k-1)*p3);
            xyz(3,k + 3) = z0 + sqrt(2/3) - sin(beta);
        end
        xyz(:,7) = [x0, y0, z0 + sqrt(2/3)]';

    case 4 % Bicapped trigonal bipyramid
        p3 = 2*pi/3;
        pp = p3/2;
        r = 1/sqrt(3);
        beta = 0.5*pi -asin(1/3) - acos(1/sqrt(3));
        r1 = cos(beta);
        p0 = 1.5*pi - pp;
        x0 = sin(pp);
        y0 = cos(pp);
        z0 = 0;
        for k = 1 : 3
            xyz(1,k) = x0 + r*cos(p0 + (k-1)*p3);
            xyz(2,k) = y0 + r*sin(p0 + (k-1)*p3);
            xyz(3,k) = z0;
        end
        xyz(:,4) = [x0 + r1*cos(p0 + pp), y0 + r1*sin(p0 + pp), z0 + sqrt(2/3) - sin(beta)]';
        xyz(:,5) = [x0 + r1*cos(p0 + pp + p3), y0 + r1*sin(p0 + pp+p3), z0 - sqrt(2/3) + sin(beta)]';
        xyz(:,6) = [x0, y0, z0 - sqrt(2/3)]';
        xyz(:,7) = [x0, y0, z0 + sqrt(2/3)]';

    otherwise % random configuration
        hR = 0.01;
        xyz = zeros(3,Na);
        xyz(:,1) = [0;0;0];
        a = randn(3,Na - 1);
        rad = sqrt(a(1,:).^2 + a(2,:).^2 + a(3,:).^2);
        a = a*diag(1./rad);
        for i = 2 : Na
            clear rad
            clear x
            rad = sqrt(xyz(1,1 : i - 1).^2 + xyz(2,1 : i - 1).^2 + xyz(3,1 : i - 1).^2);
            R = max(rad) + rstar;
            xa = R*a(:,i-1);
            x = [xyz(:,1 : i - 1), xa];
            f = LJ(x(:));
            fnew = f;
            while 1
                R = R - hR;
                xa = R*a(:,i - 1);
                x = [xyz(:,1 : i - 1), xa];
                f = fnew;
                fnew = LJ(x(:));
                if fnew > f
                    break;
                end
            end
            xyz(:,i) = xa;
        end
        cmass = mean(xyz,2);
        xyz = xyz - cmass*ones(1,Na);
end
xyz = xyz*rstar;
end

%%
%%
function v = LJpot(xyz) % 
% xyz must be a column vector
[m, ~] = size(xyz);
Na = (m + 6)/3 ;
x_aux = [0;0;0;xyz(1);0;0;xyz(2:3);0;xyz(4:m)];
m = length(x_aux);
% restore atomic coordinates 
x = reshape(x_aux,3,Na);
r = zeros(Na);
for k = 1 : Na
    r(k,:) = sqrt(sum((x - x(:,k)*ones(1,Na)).^2,1));
end
r = r + diag(ones(Na,1));
aux = 1./r.^6;
L = (aux - 1).*aux;
L = L - diag(diag(L));
v = 2*sum(sum(L));
end
 
 %%
function dv = LJgrad(xyz) % 
[m, ~] = size(xyz);
Na = (m + 6)/3 ;
x_aux = [0;0;0;xyz(1);0;0;xyz(2:3);0;xyz(4:m)];
m = length(x_aux);
% restore atomic coordinates
g = zeros(size(x_aux));
x = reshape(x_aux,3,Na);
r = zeros(Na);
for k = 1 : Na
    r(k,:) = sqrt(sum((x - x(:,k)*ones(1,Na)).^2,1));
end
r = r + diag(ones(Na,1));
L = -6*(2./r.^6 - 1)./r.^8;
for k = 1 : Na
    Lk = L(:,k);
    g(1 + (k-1)*3 : k*3) = (x(:,k)*ones(1,Na) - x)*Lk;
end

g = 4*g;
dv = g([4,7,8,10 : Na*3]);
end

%%
 function dv = LJder(r)
 dv = -12*r.^(-13) + 6*r.^(-7);
 end
%%
%%
function H = LJhess(x) % find the Hessian using finite differences
h = 1e-6;
n = length(x);
H = zeros(n);
e = eye(n);
for i = 1 : n
    di = e(:,i)*h;
    Hei = 0.5*(LJgrad(x + di) - LJgrad(x - di))/h;
    for j = 1 : i
        H(j,i) = e(j,:)*Hei;
        H(i,j) = H(j,i);
    end
end
end

%%
function x = remove_rotations_translations(xyz)
% removes rotational and translational degrees of freedom
% input should be 3 by Na matrix
% output = column vector 3*Na - 6 by 1
[m, Na] = size(xyz);
if m ~= 3 || Na < 3
    fprintf('Error in remove_rotations_translations: [m = %d, Na = %d]\n',m,Na);
    return
end
    % shift atom 1 to the origin;
    xyz = xyz - xyz(:,1)*ones(1,Na);
    % Use a Householder reflection to place atom 2 on the x-axis
    u = xyz(:,2);
    noru = norm(u);
    ind = find(abs(u) > 1e-12, 1 );
    xyz = circshift(xyz,[1 - ind,0]);
    u = xyz(:,2);
    u(1) = u(1) + sign(u(1))*noru;
    noru = norm(u);
    if noru > 1e-12
        House = eye(3);
        u = u/noru;
        House = House - 2*u*u';
        xyz = House*xyz;
    end
    % Perform rotation around the x-axis to place atom 3 onto the xy-plane
    R = eye(3);
    a = xyz(:,3);
    r = sqrt(a(2)^2 + a(3)^2);
    if r > 1e-12
        R(2,2) = a(2)/r;
        R(3,3) = R(2,2);
        R(2,3) = a(3)/r;
        R(3,2) = -R(2,3);
        xyz = R*xyz;
    end
    %% prepare input vector
    x_aux = xyz(:); % make xyz into column vector
    % x = [x2 x3 y3 x4 y4 z4 ... xNa yNa ZNa]';
    x = x_aux([4,7,8,10 : Na*3]);
%     fprintf('Ndim = %d\n',N);

end
%%
function v = LJ(xyz) % 
global feval
feval = feval + 1;
m = length(xyz);
Na = m/3;
x = reshape(xyz,3,Na);
r = zeros(Na);
for k = 1 : Na
    r(k,:) = sqrt(sum((x - x(:,k)*ones(1,Na)).^2,1));
end
r = r + diag(ones(Na,1));
aux = 1./r.^6;
L = (aux - 1).*aux;
L = L - diag(diag(L));
v = 2*sum(sum(L));
end

%% 
function drawconf(xyz,fig)
figure(fig);   
clf;
hold on;
Na = size(xyz,2);

rstar = 0.5*2^(1/6);
t = linspace(-1,1,50);
[X, Y, Z] = meshgrid(t,t,t);
V = sqrt(X.^2 + Y.^2 + Z.^2);
[faces,verts] = isosurface(X,Y,Z,V,rstar);
nverts = size(verts,1);
e = ones(nverts,1);
col = jet(Na);


hold on;
for i = 1 : Na
    patch('Vertices',verts + e*xyz(:,i)','Faces',faces,'Facecolor',col(i,:),...
        'EdgeColor','none','LineStyle','none');
end
camlight('headlight');
lighting gouraud
set(gca,'DataAspectRatio',[1,1,1]);
alpha(0.8);
view(3)
axis off
drawnow
end
