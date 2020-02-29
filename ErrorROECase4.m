%%
clc
clear all
close all

% Initial data set
data = 1;
u=0.25;
g=1;
switch data
    case 1
        hIC=@(x) 1-0.2*sin(2*pi*x);
        mIC=@(x) 0.*x;
        S=@(x,t) [0.*x;0.*x];
        bc = 'Periodic';
end

%Exact solution
load('solCase3.mat')
U_exact=U;
dx_ex=0.0002;
x_exact=0+dx_ex/2:dx_ex:2-dx_ex/2;

%definition of the number of grid point and the error vectors
M=[100,200,500,1000];
error_h=zeros(1,4);
error_m=zeros(1,4);

%compute the errors
for i=1:4
    % Define Discretization and time parameters
    dx=2.0/M(i);
    xf = 0:dx:2;
    xc = (0.5*dx):dx:(2-0.5*dx);
    N  = length(xc);
    Tfinal = 2;
    CFL    = 0.5;


    % Averaging initial conditions
    % Cell-center values sufficient for first-order schemes
    U = [hIC(xc);mIC(xc)];

    time = 0; iter = 0;

    h_ex=@(x,t) hIC(x-t);
    m_ex=@(x,t) 0.25*h_ex(x,t); 
    % Solve
    while time < Tfinal
    
        k = CFL*dx/max(abs(U(2,:)./U(1,:)) + sqrt(g*U(1,:)));
        if(time + k > Tfinal)
            k = Tfinal - time;
        end
    
        % Applying boundary conditions to obtain extended vector
        U_ext = apply_bc(U,bc);
    
        Flux = RoeFlux(U_ext(:,1:end-1),U_ext(:,2:end),g);
        U = U - k/dx*(Flux(:,2:end) - Flux(:,1:end-1));
    
        time = time + k;
        iter = iter + 1;
    end
    error_h(i)=max(abs(U(1,:)-interp1(x_exact,U_exact(1,:),xc,'linear')));
    error_m(i)=max(abs(U(2,:)-interp1(x_exact,U_exact(2,:),xc,'linear')));
end
figure(1)
loglog(1./M,error_h,'r')
hold on;
loglog(1./M,10./(M),'--k')
legend('Error h as function of dx at time 2','slope dx')

figure(2)
loglog(1./M,error_m,'r');
hold on;
loglog(1./M,10./M,'--k');
legend('Error m as function of dx at time 2','slope dx')