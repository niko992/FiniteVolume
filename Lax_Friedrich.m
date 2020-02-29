%%
clc
clear all
close all

% Initial data set
data = 4;
u=0.25;
g=1;
switch data
    case 1
        u=0.25;
        hIC =@(x) 1+0.5*sin(pi*x);
        mIC =@(x) u*hIC(x);
        S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
        bc  = 'Periodic';
    case 2
        u=0.25;
        hIC =@(x) 1+0.5*sin(pi*x);
        mIC =@(x) u*hIC(x);
        S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
        bc  = 'Open';
    case 3
        u=0.0;
        hIC=@(x) 1-0.1*sin(pi*x);
        mIC=@(x) 0.*x;
        S=@(x,t) [0.*x;0.*x];
        bc = 'Periodic';
    case 4
        hIC=@(x) 1-0.2*sin(2*pi*x);
        mIC=@(x) 0.5+0.*x;
        S=@(x,t) [0.*x;0.*x];
        bc = 'Periodic';
    case 5
        %Useful to compute convergence error of Roe method
        %Remember to change Tfinal=0.5
        hIC=@(x) 1+0.*x;
        mIC=@(x) -1.5.*(x<=1)+0.*(x>1);
        S=@(x,t) [0.*x;0.*x];
        bc = 'Open';
end


dx=0.02;
xf = 0:dx:2;
xc = (0.5*dx):dx:(2-0.5*dx);
N  = length(xc);
Tfinal = 2.;
CFL    = 0.5;

if data~=1
    myfilename=strcat('solCase',num2str(data),'.mat');
    load(myfilename);
    U_exact=U;
    x_exact=0.0001:0.0002:1.9999;
end
% Averaging initial conditions
% Cell-center values sufficient for first-order schemes
U = [hIC(xc);mIC(xc)];

time = 0; iter = 0;

%source term and exact solutions
%S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
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
    
    Flux = LaxFriedFlux(U_ext(:,1:end-1),U_ext(:,2:end),k,g,dx);
    U = U - k/dx*(Flux(:,2:end) - Flux(:,1:end-1)) + k*(S(xc,time));
    
    time = time + k;
    iter = iter + 1;
    
        
end
if data == 1
    figure(1)
    subplot(2,1,1)
    plot(xc,U(1,:),'-r','LineWidth',2);
    hold on;
    plot(xc,h_ex(xc,time),'--k');
    legend('h obtained by LF','true solution h');
    hold off;
      
    title(['Solution with dx=',num2str(dx), ' and Time = ',num2str(time)])

    subplot(2,1,2)
    plot(xc,U(2,:),'-r','LineWidth',2);
    hold on;
    plot(xc,m_ex(xc,time),'--k');
    legend('m obtained by LF','true solution m');
    hold off;
    ylim([-0.5 2.3]);xlim([0 2]);
else
    figure(1)
    subplot(2,1,1)
    plot(xc,U(1,:),'-r','LineWidth',2);
    hold on
    plot(xc,interp1(x_exact,U_exact(1,:),xc,'linear'))
    legend('h obtained by LF','exact solution');
      
    title(['Solution with dx=',num2str(dx), ' and Time = ',num2str(time)])
    hold off
    
    subplot(2,1,2)
    plot(xc,U(2,:),'-r','LineWidth',2);
    hold on
    plot(xc,interp1(x_exact,U_exact(2,:),xc,'linear'))
    legend('m obtained by LF','exact solution');
    hold off;
    ylim([-1.1 2]);xlim([0 2]);    
end


