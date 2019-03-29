%% What we are doing???
% Volpi rabbiose
% 2 dimension
% In this case we choose Neumann BCs.
%% SPACE DOMAIN
L1 = 100;
L2 = 1;
Deltax = 0.1;
Deltay = 0.1;
Nx = L1/Deltax-1;
Ny = L2/Deltay-1;
s1 = linspace(0,L1,Nx+2);
s2 = linspace(0,L2,Ny+2);
s1 = s1(2:Nx+1);
s2 = s2(2:Ny+1);

%% TIME DOMAIN
Deltat = 0.2;
Tmax = 180;
Tout = 1;

%% WRITE THE APPROXIMATE SCHEME (Crank-Nicholson method)
d = 1.0; % from 60 to 330 km^2 per year
alpha = 0.5*Deltat/Deltax^2*d*d;
beta = 0.5*Deltat/Deltay^2*d*d;
% generate matrix of tridiagonal system
aux = zeros(Ny*Nx,5);
aux(:,3) = 1 + 2*alpha + 2*beta  ; % main diagonal

for j = 1:Ny
    aux((j-1)*Nx+1,3) = aux((j-1)*Nx+1,3) - beta;
    aux(j*Nx,3) = aux(j*Nx,3) - beta;
end
for i = 1:Nx
    aux(i,3) = aux(i,3) - alpha;
    aux((Ny-1)*Nx+i,3) = aux((Ny-1)*Nx+i,3) - alpha;
end

aux(:,2) = -alpha;
aux(:,1) = -beta;
aux(:,4) = -alpha;
aux(:,5) = -beta;
for j = 1:Ny
    aux(j*Nx,2) = 0;
    aux((j-1)*Nx+1,4) = 0;
end
A = (spdiags( aux, [-Nx -1 0 1 Nx], Ny*Nx,Ny*Nx))'; % transposed of the tridiagonal matrix with the given values


aux(:,3) =  1 - 2*alpha - 2*beta ; % main diagonal
for j = 1:Ny
    aux((j-1)*Nx+1,3) = aux((j-1)*Nx+1,3) + beta;
    aux(j*Nx,3) = aux(j*Nx,3) + beta;
end
for i = 1:Nx
    aux(i,3) = aux(i,3) + alpha;
    aux((Ny-1)*Nx+i,3) = aux((Ny-1)*Nx+i,3) + alpha;
end
aux(:,2) = alpha;
aux(:,1) = beta;
aux(:,4) = alpha;
aux(:,5) = beta;
for j = 1:Ny
    aux(j*Nx,2) = 0;
    aux((j-1)*Nx+1,4) = 0;
end
B = (spdiags( aux, [-Nx -1 0 1 Nx], Ny*Nx,Ny*Nx))'; % transposed of the tridiagonal matrix with the given values

%% PARAMETERS
a = 1.0; %1 per year
b = 0.5; %0.5 per year
k = 1; %0.25 to 4 foxes per km^-2
bet = 80.0; % 80 km^2 per year
sigma = 1.0/28; % days^-1
alfa = 1.0/5; % days^-1
cmat = zeros(Nx, Ny);
cmat(1:floor(Nx/4),:) = a*7/3;
c = reshape(cmat,Ny*Nx,1);
%Kt = (sigma + a)*(alfa + a)/(sigma*bet);

%% 4 SPECIES
S = zeros(Ny*Nx, 1);
S(:) = k;
I = zeros(Ny*Nx, 1);
Rmat = zeros(Nx, Ny);
Rmat(Nx-1,:) = 0.001;
R = reshape(Rmat,Ny*Nx,1);
V = zeros(Ny*Nx, 1);

M = Tmax/Deltat;
memory = cell(3,M+1);
t = 0;
mat = reshape(R, [Nx, Ny]);

memory{1,1} = mat(:,floor(Ny/2));
[L,U,P] = lu(A);
for j=1:M % temporal steps
    t = t + Deltat;
    N = S + I + R + V;
    SN = S + ((-c+(a-b)*(1-N/k) - bet*R).*S+a*V)*Deltat;
    IN = I + (bet*S.*R - (sigma + b +(a-b)*N/k).*I)*Deltat;
    rhs = B*R + (sigma*I - (alfa + b + (a-b)*N/k).*R)*Deltat;
    y = L\rhs;
    RN = U\y; % solution through two triangular systems
    VN = V + (c.*S-(b+(a-b)/k*N).*V)*Deltat;
    S = SN;
    I = IN;
    R = RN;
    V = VN;
    mat = reshape(R, [Nx, Ny]);
    memory{1,j} = flipud(mat(:,floor(Ny/2)));
end

%% 3 SPECIES K
kmat = ones(Nx, Ny);
kmat(1:floor(Nx/4),:) = 0.3;
k = reshape(kmat,Ny*Nx,1);

S = k;
I = zeros(Ny*Nx, 1);
Rmat = zeros(Nx, Ny);
Rmat(Nx-1,:) = 0.001;
R = reshape(Rmat,Ny*Nx,1);
V = zeros(Ny*Nx, 1);

M = Tmax/Deltat;
t = 0;
mat = reshape(R, [Nx, Ny]);
memory{2,1} = mat(:,floor(Ny/2));
for j=1:M % temporal steps
    t = t + Deltat;
    N = S + I + R;
    SN = S + ((a-b)*(1-N./k) - bet*R).*S*Deltat;
    IN = I + (bet*S.*R - (sigma + b +(a-b)*N./k).*I)*Deltat;
    rhs = B*R + (sigma*I - (alfa + b + (a-b)*N./k).*R)*Deltat;
    y = L\rhs;
    RN = U\y; % solution through two triangular systems
    S = SN;
    I = IN;
    R = RN;
    mat = reshape(R, [Nx, Ny]);
    memory{2,j} = flipud(mat(:,floor(Ny/2)));
end

%% 3 SPECIES BETA
k = 1.0;
betmat = ones(Nx, Ny)*80.0; % 80 km^2 per year
betmat(1:floor(Nx/4),:) = 80.0*0.3;
bet = reshape(betmat,Ny*Nx,1);

S = zeros(Ny*Nx, 1);
S(:) = k;
I = zeros(Ny*Nx, 1);
Rmat = zeros(Nx, Ny);
Rmat(Nx-1,:) = 0.001;
R = reshape(Rmat,Ny*Nx,1);
V = zeros(Ny*Nx, 1);

M = Tmax/Deltat;
t = 0;
mat = reshape(R, [Nx, Ny]);
memory{2,1} = mat(:,floor(Ny/2));
for j=1:M % temporal steps
    t = t + Deltat;
    N = S + I + R;
    SN = S + ((a-b)*(1-N/k) - bet.*R).*S*Deltat;
    IN = I + (bet.*S.*R - (sigma + b +(a-b)*N/k).*I)*Deltat;
    rhs = B*R + (sigma*I - (alfa + b + (a-b)*N/k).*R)*Deltat;
    y = L\rhs;
    RN = U\y; % solution through two triangular systems
    S = SN;
    I = IN;
    R = RN;
    mat = reshape(R, [Nx, Ny]);
    memory{3,j} = flipud(mat(:,floor(Ny/2)));
end

%% PLOT BARRIERA
fig = figure;
i=1;
t=0;
for j=1:M
    t = t + Deltat;
    if (t >= i*Tout)
        i = i + 1;
        p1 = plot(memory{1,j});
        hold on;
        p2 = plot(memory{2,j});
        p3 = plot(memory{3,j});
        hold off;
        ylim([0, 0.005]);
        title('4');
        h = [p1;p2;p3];
        legend(h,'4 specie','barriera k','barriera beta');
    end
end

%% COSE
barriere = zeros(3, Nx);
for j=1:M
    for i=1:3
        for x=1:Nx
            barriere(i,x) = max(barriere(i,x), memory{i,j}(x));
        end
    end
end

fig = figure;
p1 = plot(barriere(1,:));
hold on;
p2 = plot(barriere(2,:));
p3 = plot(barriere(3,:));
hold off;
xlim([145*5, 165*5]);
ylim([0, 0.008]);
title('barriere');
h = [p1;p2;p3];
legend(h,'4 specie','barriera k','barriera beta');