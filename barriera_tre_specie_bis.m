%% What we are doing???
% Volpi rabbiose
% 2 dimension
% In this case we choose Neumann BCs.
%% SPACE DOMAIN
L1 = 4;
L2 = 140;
Deltax = 0.2;
Deltay = 0.2;
Nx = L1/Deltax-1;
Ny = L2/Deltay-1;
s1 = linspace(0,L1,Nx+2);
s2 = linspace(0,L2,Ny+2);
s1 = s1(2:Nx+1);
s2 = s2(2:Ny+1);

%% TIME DOMAIN
Deltat = 0.02;
Tmax = 300;
Tout = 3;

%% PARAMETERS
a = 1.0; %1 per year
b = 0.5; %0.5 per year
kmat = ones(Nx, Ny)*2;
kmat(:,floor(80/Deltax):Ny) = 0.4;
k = reshape(kmat,Ny*Nx,1);
bet = 80.0; % 80 km^2 per year
sigma = 1.0/28; % days^-1
alfa = 1.0/5; % days^-1
d = 1.0; % from 60 to 330 km^2 per year

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
B = (spdiags(aux, [-Nx -1 0 1 Nx], Ny*Nx,Ny*Nx))'; % transposed of the tridiagonal matrix with the given values

%% Initial Conditions
S = k;
I = zeros(Ny*Nx, 1);
Rmat = zeros(Nx, Ny);
Rmat(:,2) = 0.001;
R = reshape(Rmat,Ny*Nx,1);

i = floor(Tmax/Tout);
M = Tmax/Deltat;
memory = cell(3,i+1);
t = 0;
temp = reshape(S, [Nx, Ny]);
memory{1,1} = temp(floor(Nx/2),:);
temp = reshape(I, [Nx, Ny]);
memory{2,1} = temp(floor(Nx/2),:);
temp = reshape(R, [Nx, Ny]);
memory{3,1} = temp(floor(Nx/2),:);

contagiati_c = zeros(Ny, 1);
contagiati_t = zeros(Ny, 1);
memory{4,1} = contagiati_c;
integrale_c = zeros(i+1);

%% TIME EVOLUTIONARY STEPS
[L,U,P] = lu(A);
i=1;
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
    for m=1:Ny
        index = floor(m*Nx-Nx/2);
        if (R(index)+I(index)>contagiati_c(m))
            contagiati_c(m) = R(index)+I(index);
            contagiati_t(m) = t;
        end
    end
    if (t >= i*Tout)
        i = i + 1;
        temp = reshape(S, [Nx, Ny]);
        memory{1,i} = temp(floor(Nx/2),:);
        temp = reshape(I, [Nx, Ny]);
        memory{2,i} = temp(floor(Nx/2),:);
        temp = reshape(R, [Nx, Ny]);
        memory{3,i} = temp(floor(Nx/2),:);
        memory{4,i} = contagiati_c;
        integrale_c(i) = 0;
        
        for m=floor(80/Deltax):Ny
            index = floor(m*Nx-Nx/2);
            integrale_c(i) = integrale_c(i) + (R(index)+I(index))/L2;
        end
    end
end

%% PLOT SIMULAZIONE
fig = figure;
pause
for j=1:i
    subplot(1,2,1);
    plot(1:Ny, memory{1,j});
    ylim([0, 1]);
    title('Volpi sane');
    
    subplot(1,2,2);
    plot(1:Ny, memory{2,j}+memory{3,j});
    ylim([0, 1]);
    hold on;
    plot(1:Ny, memory{4,j});
    ylim([0, 1]);
    title('Volpi infette');
    hold off;
    set(fig, 'NumberTitle', 'off', 'Name', sprintf('Time: %f', (j-1)*Tout));
    pause(0.05)
end

%% PLOT BARRIERA
fig = figure;
plot(integrale_c);
title('integrale nel tempo')

%% GIF
fig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for j=1:i
    % Draw plot for y = x.^n
    subplot(1,2,1);
    plot(1:Ny, memory{1,j});
    ylim([0, 1]);
    title('Volpi sane');
    
    subplot(1,2,2);
    plot(1:Ny, memory{2,j}+memory{3,j});
    ylim([0, 1]);
    hold on;
    plot(1:Ny, memory{4,j});
    ylim([0, 1]);
    title('Volpi infette');
    hold off;
    set(fig, 'NumberTitle', 'off', 'Name', sprintf('Time: %f', (j-1)*Tout));
    drawnow 

    % Capture the plot as an image 
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 

    % Write to the GIF File 
    if j == 1
      imwrite(imind,cm,filename,'gif', 'DelayTime',0.1, 'Loopcount',inf);
    else
      imwrite(imind,cm,filename,'gif', 'DelayTime',0.1,'WriteMode','append');
    end
end