%% What we are doing???
% Volpi rabbiose
% 2 dimension
% In this case we choose Neumann BCs.
%% SPACE DOMAIN
L1 = 60;
L2 = 60;
Deltax = 0.5;
Deltay = 0.5;
Nx = L1/Deltax-1;
Ny = L2/Deltay-1;
s1 = linspace(0,L1,Nx+2);
s2 = linspace(0,L2,Ny+2);
s1 = s1(2:Nx+1);
s2 = s2(2:Ny+1);

%% TIME DOMAIN
Deltat = 0.01;
Tmax = 300;
Tout = 3;

%% WRITE THE APPROXIMATE SCHEME (Crank-Nicholson method)
alpha = 0.5*Deltat/Deltax^2;
beta = 0.5*Deltat/Deltay^2;
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
%0.25 to 4 foxes per km^-2
kmat = ones(Nx, Ny);
kmat(50:70,25:35) = 2;
kmat(10:35,:) = 0.2;
k = reshape(kmat,Ny*Nx,1);
%beta = 80.0; % 80 km^2 per year
betamat = zeros(Nx, Ny);
betamat(:,:) = 80.0;
betamat(:,85:110) = 8.0;
beta = reshape(betamat,Ny*Nx,1);
sigma = 1.0/28; % days^-1
alfa = 1.0/5; % days^-1
d = 100.0; % from 60 to 330 km^2 per year

%% Initial Conditions
S = k;
I = zeros(Ny*Nx, 1);
Rmat = zeros(Nx, Ny);
Rmat(Nx-1,2) = 0.001;
R = reshape(Rmat,Ny*Nx,1);

M = Tmax/Deltat;
memory = cell(3,M+1);
t = 0;
memory{1,1} = reshape(S, [Nx, Ny]);
memory{2,1} = reshape(I, [Nx, Ny]);
memory{3,1} = reshape(R, [Nx, Ny]);

%% TIME EVOLUTIONARY STEPS
[L,U,P] = lu(A);
i=1;
for j=1:M % temporal steps
    t = t + Deltat;
    N = S + I + R;
    SN = S + ((a-b)*(1-N./k) - beta.*R).*S*Deltat;
    IN = I + (beta.*S.*R - (sigma + b +(a-b)*N./k).*I)*Deltat;
    rhs = B*R + (sigma*I - (alfa + b + (a-b)*N./k).*R)*Deltat;
    y = L\rhs;
    RN = U\y; % solution through two triangular systems
    S = SN;
    I = IN;
    R = RN;
    if (t >= i*Tout)
        i = i + 1;
        memory{1,i} = reshape(S, [Nx, Ny]);
        memory{2,i} = reshape(I, [Nx, Ny]);
        memory{3,i} = reshape(R, [Nx, Ny]);
    end
end

%% PLOT
fig = figure;
pause
for j=1:i
    subplot(1,3,1)
    surfl(s2, s1, memory{1,j});
    zlim([0, 2])
    colormap(pink);
    shading interp;
    title('Volpi sane')
    
    subplot(1,3,2)
    surfl(s2, s1, memory{2,j});
    zlim([0, 2])
    colormap(pink);
    shading interp;
    title('Volpi quiescenti')
    
    subplot(1,3,3)
    surfl(s2, s1, memory{3,j});
    zlim([0, inf])
    colormap(pink);
    shading interp;
    title('Volpi rabbiose')
    
    set(fig, 'NumberTitle', 'off', 'Name', sprintf('Time: %f', (j-1)*Tout));
    pause(0.05)
end

