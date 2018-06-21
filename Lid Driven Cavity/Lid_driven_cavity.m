% Lid Driven Cavity Problem [Validated with Ghia et al]
% @author: YOGESHWARAN R

clear;
clc;

%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRID GENERATION
l=1.0; h=1.0;       % Square duct
nx=31;ny=31;        % No. of Grids

j1=20; j2=25;       % Left Opening
j3=5;j4=10;         % Right Opening

% TIME STEPPING
dt=0.001; 

% OTHER IMP. COMPUTATIONAL SPECIFICATIONS
up = 1.0;                     % Free Stream Velocity on top
Re = 1000;                    % Renold's Number
rf = 1.5;                     % Relaxation Factor
sf1 = 0.1;                    % arbitrary constant
maxerr=0.001;                   
MaxIter=1000;
err1 = 1000.0;
itc = 0;

%%%%%%%%%%%%%% DOMAIN DISCRETIZATION AND INITIALIZATION %%%%%%%%%%%%%%%
dx=l/(nx-1);
dy=h/(ny-1);

%Initialize the Parameters
sf=zeros(nx,ny);
sfold=zeros(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);
w=zeros(nx,ny);
rw=zeros(nx,ny);
rsf=zeros(nx,ny);

for i=1:nx
    u(i,ny) = up;       % Top wall moving with velocity "up"
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(err1>0.002)
       
    wold = w;
    %%% SF B.C %%%%%
    for i=1:nx
        sf(i,ny) = 0;
    end       
        
    for iter=1:1000                
        %SOLVING CONTINUITY EQUATION FOR OBTAINING STREAM FUNCTION
        sfold = sf;
        for i=2:nx-1                        
            for j=2:ny-1
                sf(i,j) = (rf/4.0)*((w(i,j)*dx*dx)+sf(i+1,j)+ ...
                    sf(i-1,j)+sf(i,j+1)+sf(i,j-1))+((1.-rf) *sf(i,j));
            end
        end
        err=0.0;
        for i =1:nx;
            for j = 1:ny,
                err =  err + abs(sfold(i,j)-sf(i,j));
            end;
        end;
        if err<=0.001            
            break;
        end
    end
    
    %%%%%%%%%%%%  VORTICITY BOUNDARY CONDITION %%%%%%%%%%%%%%%%%  
    for i=1:nx
        %Top B.C
        w(i,ny) = -2.0*((sf(i,ny-1)-sf(i,ny))/(dy*dy) + (up/dy));
        
        %Bottom B.C
        w(i,1) = 2* (sf(i,1)  - sf(i,2)  ) / (dy*dy);
    end
    
    % Left B.C
    for j=2:j1
        w(1,j) = 2.0 * (sf(1,j) - sf(2,j) ) / (dx*dx);
    end
    
    for j=j1+1:j2
        w(1,j) = 2.0 * (sf(1,j) - sf(2,j) ) / (dx*dx);
    end
    
    for j=j2+1:ny-1
        w(1,j) = 2.0 * (sf(1,j) - sf(2,j) ) / (dx*dx);
    end
        
    % Right B.C
    for j =2:j3
        w(nx,j)= 2.0 * ((sf(nx,j)-sf(nx-1,j)) / (dx*dx));
    end
    
    for j =j3+1:j4
        w(nx,j)= 2.0 * ((sf(nx,j)-sf(nx-1,j)) / (dx*dx));
    end
        
    for j =j4+1:ny-1
        w(nx,j)= 2.0 * ((sf(nx,j)-sf(nx-1,j)) / (dx*dx));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALC U AND V
    for i=2:nx-1
        for j=2:ny-1
            u(i,j) = 0.5 * (sf(i,j+1) - sf(i,j-1) ) / (dy);
            v(i,j) = -0.5 *(sf(i+1,j) - sf(i-1,j) )/ (dx);
        end
    end    
   
    % SOLVING VORTICITY EQUATION
    for i=2:nx-1
        for j=2:ny-1
            term1 = (w(i+1,j)-2.0*w(i,j)+w(i-1,j) ) / (dx*dx);
            term2 = (w(i,j+1)-2.0*w(i,j)+w(i,j-1) ) / (dy*dy);
            term3 = (u(i+1,j)*w(i+1,j) - u(i-1,j)*w(i-1,j) )/ (2.0*dx);
            term4 = (v(i,j+1)*w(i,j+1) - v(i,j-1)*w(i,j-1) )/ (2.0*dy);            
            rw(i,j) =   (term1+term2)/Re-term3-term4;
            w(i,j) = w(i,j) + dt*rw(i,j);
        end
    end
    
    err1=0.0;
    for i =1:nx;
        for j = 1:ny,
            err1 =  err1 + abs(wold(i,j)-w(i,j));
        end;
    end;   
    
    itc = itc+1    
end
%%
figure(1)
subplot(1,2,1)
contour(w','DisplayName','U1','showtext','on','LineWidth',2.0,'LevelStep',1.0);  % adjust levelstep for diff. Re
title('Vorticity')
xlabel('X') 
ylabel('Y')
colorbar('southoutside')

subplot(1,2,2)
contour(sf','DisplayName','U1','showtext','on','LineWidth',1.8,'LevelStep',0.005);
title('Stream Function')
xlabel('X') 
ylabel('Y')
colorbar('southoutside')

figure(2)
contour(u','DisplayName','U1','showtext','on','LineWidth',2.0,'LevelStep',0.1);
title('U - Contour')
xlabel('X') 
ylabel('Y')
colorbar
hold on
quiver(u',v',1.5,'--k')
hold off

