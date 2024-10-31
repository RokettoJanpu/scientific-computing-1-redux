function [A,b,sol] = TestMatrixA()
flag = 1;
nx = 5;
ny = 6;
% Solving Equation -grad * (a(x,y) * grad(u)) = f(x,y)
% p = 9;
% n = 2^p + 1; % the finest mesh is n-by-n
% h = 2^(-p);
nx1 = nx - 1;
ny1 = ny - 1;
hx = 1/nx1;
hy = 1/ny1;
Ix = 2 : nx1;
Iy = 2 : ny1;

tx = linspace(0,1,nx);
ty = linspace(0,1,ny);
[x,y] = meshgrid(tx,ty);
% set f(x,y) so that the exact solution is the given one 
uexact = x.*(1-x).*y.*(1-y); % u exact
a = 1 + x + 2*y.^2;
f = diff_oper(hx,hy,uexact,a);

A = make_matrix(a,hx,hy);
u = zeros(ny,nx);  % the initial guess for the solution
b = f(:);

if flag == 1
    u(Iy,Ix) = reshape(A\b,ny-2,nx-2);
    fprintf('Max Error = %d\n',max(max(abs(u - uexact))));
end

if flag == 2 % use preconditioned CG
    u(Iy,Ix) = reshape(CG(A,b),ny-2,nx-2);
end

uinner = u(Iy,Ix);
sol = uinner(:);
end
%%
function Lu = diff_oper(hx,hy,u,a)
[ny,nx] = size(u);

    as = 0.5*(a + circshift(a,[1,0]))/hy^2;
    an = 0.5*(a + circshift(a,[-1,0]))/hy^2;
    aw = 0.5*(a + circshift(a,[0,1]))/hx^2;
    ae = 0.5*(a + circshift(a,[0,-1]))/hx^2;
    ap = aw + ae + as + an;
    Ix = 2 : nx - 1;
    Iy = 2 : ny - 1;

Lu = ap(Iy,Ix).*u(Iy,Ix) - as(Iy,Ix).*u(Iy - 1,Ix) ...
    - an(Iy,Ix).*u(Iy + 1,Ix) - aw(Iy,Ix).*u(Iy,Ix - 1)...
    - ae(Iy,Ix).*u(Iy,Ix + 1);

end

%%
function A = make_matrix(a,hx,hy)
hx2 = hx^2;
hy2 = hy^2;
[ny,nx] = size(a);
ny2 = ny - 2;
nx2 = nx - 2;
cN = 0.5*(a + circshift(a,[-1 0]))/hy2;
cS = 0.5*(a + circshift(a,[1 0]))/hy2;
cW = 0.5*(a + circshift(a,[0 1]))/hx2;
cE = 0.5*(a + circshift(a,[0 -1]))/hx2;
cP = -(cE + cW + cN + cS);
    Ix = 2 : nx - 1;
    Iy = 2 : ny - 1;
an = -cN(Iy,Ix);
as = -cS(Iy,Ix);
aw = -cW(Iy,Ix);
ae = -cE(Iy,Ix);
ap = -(an + as + ae + aw);
% Conversion to column vectors
cn = an(:); cs = as(:); cw = aw(:); ce = ae(:); cp = ap(:);

% First set up the mathrix A without taking care of B and D
% A has nx2-by-nx2 blocks of size ny2-by-ny2
A = spdiags([circshift(cw,[-ny2 0]),circshift(cs,[-1 0]),cp,...
    circshift(cn,[1 0]),circshift(ce,[ny2 0])],[-ny2,-1,0,1,ny2],ny2*nx2,ny2*nx2);
% Take care of the Dirichlet BCs: u_y(x,0) = u_y(x,1) = 0
for j = 1 : nx2 - 1 
    A(ny2*j,ny2*j + 1) = 0;
    A(ny2*j + 1,ny2*j) = 0;
end
end
%%

function x = CG(A,b)
%% solve Ax = b using CG
n = length(A); % the size of the matrix A

figure(1); clf; hold on;
col = ['r','b'];
set(gca,'YScale','log','Fontsize',20);
xlabel('Iter #','Fontsize',20);
ylabel('Residual','Fontsize',20);

x0 = zeros(length(A),1);
%%
for run = 1 : 2
    if run == 1
        % np preconditioning
        M = speye(n);
        L = speye(n);
    else       
        % Incomplete Cholesky preconditioning
        L = ichol(A);
        M = L*L';
    end
    clear nr
    x = x0;
    tol = 1e-9;
    r = A*x-b;
    z = L\r;
    y = L'\z;
%     y = M\r;
    p = -y;
    k = 1;
    ry = r'*y;
    nor = norm(r);
    fprintf('k = %d\t res = %d\n',k,nor);
    nr(1) = nor;
    plot(k,nor,'.','Markersize',10,'color',col(run));
    tic
    while nor > tol
        Ap = A*p;
        al = ry/(p'*Ap);
        x = x + al*p;
        r = r + al*Ap;
        z = L\r;
        y = L'\z;
%        y = M\r;
        rynew = r'*y;
        bet = rynew/ry;
        p = -y + bet*p;
        k = k+1;
        ry = rynew;
        nor = norm(r);
        nr(k) = nor;
        if mod(k,100) == 0
            fprintf('k = %d\t res = %d\n',k,nor);
        end
    end
    plot([1:k],nr,'.','Markersize',10,'color',col(run));
    CPUtime = toc;
    fprintf('Run %d: %d iterations, CPUtime = %d\n',run,k,CPUtime);
end


c = condest((A));
fprintf('cond(A) = %d\n',c);
% figure(2); clf; hold on
% plot([1:n],sort(eig(full(A)),'descend'),'.b','MArkersize',20);
end    
  