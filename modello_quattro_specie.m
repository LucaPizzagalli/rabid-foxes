%% What we are doing???
% Volpi rabbiose
% 2 dimension
% In this case we choose Neumann BCs.
%% SPACE DOMAIN
L1 = 100;
L2 = 100;
Deltax = 0.5;
Deltay = 0.5;
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

%% PARAMETERS
a = 1.0; %1 per year
b = 0.5; %0.5 per year
k = 1; %0.25 to 4 foxes per km^-2
bet = 80.0; % 80 km^2 per year
sigma = 1.0/28; % days^-1
alfa = 1.0/5; % days^-1
d = 1.0; % from 60 to 330 km^2 per year
cmat = zeros(Nx, Ny);
cmat(1:floor(Nx/4),:) = 2.33;
c = reshape(cmat,Ny*Nx,1);
%Kt = (sigma + a)*(alfa + a)/(sigma*bet);

%% WRITE THE APPROXIMATE SCHEME (Crank-Nicholson method)
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

%% Initial Conditions
S = zeros(Ny*Nx, 1);
S(:) = k;
I = zeros(Ny*Nx, 1);
Rmat = zeros(Nx, Ny);
Rmat(Nx-1,floor(Ny/2)) = 0.001;
R = reshape(Rmat,Ny*Nx,1);
V = zeros(Ny*Nx, 1);

M = Tmax/Deltat;
memory = cell(3,M+1);
t = 0;
memory{1,1} = reshape(S, [Nx, Ny]);
memory{2,1} = reshape(I, [Nx, Ny]);
memory{3,1} = reshape(R, [Nx, Ny]);
memory{4,1} = reshape(V, [Nx, Ny]);

%% TIME EVOLUTIONARY STEPS
[L,U,P] = lu(A);
i=1;
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
    if (t >= i*Tout)
        i = i + 1;
        memory{1,i} = reshape(S, [Nx, Ny]);
        memory{2,i} = reshape(I, [Nx, Ny]);
        memory{3,i} = reshape(R, [Nx, Ny]);
        memory{4,i} = reshape(V, [Nx, Ny]);
    end
end

%% GIF
fig = figure('pos',[20 10 1200 600]);
filename = 'quattro_specie_3D.gif';
ha = tight_subplot(1,3,[.01 .04],[.1 .1],[.04 .01]);
for j=1:i
    axes(ha(1));
    surfl(s2, s1, memory{1,j});
    xlim([0, L2])
    ylim([0, L1])
    zlim([0, 2])
    %set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    %set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    shading interp;
    title('Volpi sane')
    colormap(gca, winter);
    
    axes(ha(2));
    surfl(s2, s1, memory{3,j});
    xlim([0, L2])
    ylim([0, L1])
    zlim([0, 0.04])
    %set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    %set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    shading interp;
    title('Volpi rabbiose')
    colormap(gca, spring);
    
    axes(ha(3));
    surfl(s2, s1, memory{4,j});
    xlim([0, L2])
    ylim([0, L1])
    zlim([0, 2])
    %set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    %set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    shading interp;
    title('Volpi vaccinate')
    colormap(gca, autumn);
    delete(findall(gcf,'type','annotation'));
    t = annotation('textbox', [0.5, 0.9, 0.1, 0.1], 'string', ['Time: ', num2str((j-1)*Tout)]);
    t.FontSize = 12;
    t.FontWeight = 'bold';
    set(fig, 'NumberTitle', 'off', 'Name', sprintf('Time: %f', (j-1)*Tout));
    drawnow
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    if j == 1
      imwrite(imind,cm,filename,'gif', 'DelayTime',0.05, 'Loopcount',inf);
    else
      imwrite(imind,cm,filename,'gif', 'DelayTime',0.05,'WriteMode','append');
    end
end