%%
clear all;
close all;

nen = 2;
ngp = 3;

xmin = 0;
xmax = 10000;

ne=50;
n=ne+1;
dx = (xmax-xmin)/ne;
coordx = transpose(linspace(xmin,xmax,n));

tmin = 0;
tmax = 20e13;
nt = 50;
coordt = transpose(linspace(tmin,tmax,nt));
dt = (tmax-tmin)/(nt-1);

%%

[ GQweight, GQpoint ] = gauss_quad(ngp);
GQpoint  = double(GQpoint);
GQweight = double(GQweight);

basis = zeros(nen, ngp);
dbasis = zeros(nen, ngp);

for i = 1:ngp
    ksi = GQpoint(i);
    basis(1,i)  = 0.5 * (1-ksi); % shape functions at GQ points
    basis(2,i)  = 0.5 * (1+ksi);
    dbasis(1,i) = -0.5;
    dbasis(2,i) = 0.5;
end

%% INPUT PARAMETERS

phi = 0.25;
rhor = 2650;
s = '-0.8e-6';

BCL = 10;
BCR = 225;

discont = 0;
intrusion = zeros(n,1);
intrusion(120:130,1) = 1000;


%% Residual

%Initial guess
T_final = zeros(n,nt);
T_final(:,1) = interp1([xmin, xmax],[BCL, BCR],coordx);
T0 = T_final(:,1);

iter_max=100;
norms = zeros(iter_max,nt);
maxval = 10000;
tol = 1e-8;

M1=1;


%%

fileID = fopen('fluid_heatcap.txt','r');
cf_raw = (fscanf(fileID,'%f',[2 Inf]))';
fclose(fileID);

fileID = fopen('rock_heatcap.txt','r');
cr_raw = (fscanf(fileID,'%f',[2 Inf]))';
fclose(fileID);

fileID = fopen('fluid_thermcon.txt','r');
kappaf_raw = (fscanf(fileID,'%f',[2 Inf]))';
fclose(fileID);

fileID = fopen('rock_thermcon.txt','r');
kappar_raw = (fscanf(fileID,'%f',[2 Inf]))';
fclose(fileID);

fileID = fopen('fluid_density.txt','r');
rhof_raw = (fscanf(fileID,'%f',[2 Inf]))';
fclose(fileID);


%%

tic

for l=2:nt
    
    iter=1; %reset iteration value for each time step
    T0 = T_final(:,l-1); %initial guess is the previous time step's solution
    
    %%
    
    while iter < iter_max
        
        res = zeros(n,1);
        J = zeros(n,n);
        
        for e=1:ne
            
            res_ele = zeros(nen,1);
            j_ele = zeros(nen,nen);
            
            for k = 1:ngp
                
                ksi = GQpoint(k);
                x1e = coordx(e);
                x2e = coordx(e+1);
                x = 0.5 * (x2e-x1e) * ksi + 0.5 * (x1e + x2e);
                jacob = 0.5 * (x2e-x1e);
                
                %%
                
                S = eval(s)*(1-phi);
                
                T = T0(e)*basis(1,k)+T0(e+1)*basis(2,k);
                Tdash = (T0(e)*dbasis(1,k)+T0(e+1)*dbasis(2,k))/jacob;
                T_prev = T_final(e,l-1)*basis(1,k)+T_final(e+1,l-1)*basis(2,k);
                
                F = phi*interp1(rhof_raw(:,1),rhof_raw(:,2),max(0,min(1000,T_prev)))*interp1(cf_raw(:,1),cf_raw(:,2),max(0,min(1000,T_prev))) + ...
                    (1-phi)*rhor*interp1(cr_raw(:,1),cr_raw(:,2),max(0,min(1000,T_prev)));
                
                kappaf = interp1(kappaf_raw(:,1),kappaf_raw(:,2),max(0,min(1000,T_prev)));
                kappar = interp1(kappar_raw(:,1),kappar_raw(:,2),max(0,min(1000,T_prev)));
                kappa = kappaf^phi * kappar^(1-phi);
                
                for i = 1:nen
                    
                    res_ele(i) = res_ele(i) + ( 1/dt*basis(i,k)*F*(T - T_prev)+ ...
                        kappa*Tdash*dbasis(i,k)/jacob + S*basis(i,k)) * jacob * GQweight(k);
                    
                    for j = 1:nen
                        
                        j_ele(i,j) = j_ele(i,j) + ...
                            (1/dt*basis(i,k)*basis(j,k)*F + kappa* dbasis(j,k)/jacob *dbasis(i,k)/jacob ) * jacob * GQweight(k);
                    end
                end
            end
            
            res(e)   = res(e) + res_ele(1);
            res(e+1) = res(e+1) + res_ele(2);
            J(e,e) = J(e,e) + j_ele(1,1);
            J(e,e+1) = J(e,e+1) + j_ele(1,2);
            J(e+1,e) = J(e+1,e) + j_ele(2,1);
            J(e+1,e+1) = J(e+1,e+1) + j_ele(2,2);
        end
        
        % IMPOSE BOUNDARY CONDITIONS
        res(1) = T0(1)-BCL;
        res(end) = T0(end)- BCR;
        J(1,:) = 0;
        J(1,1) = 1;
        J(n,:)  = 0;
        J(n,n) = 1;
        
        %%
        % solve for du
        dT = J\(-res);
        
        norms(iter,l) = norm(res)/sqrt(ne);
        
        if norms(iter,l) < tol
            T_final(:,l) = T0;
            
            if (all(T_final(:,l) - T_final(:,l-1) == zeros(n,1) )) && (discont == 0)
                
                T_final(26:28,l) = 1000;
                
                discont = discont + 1;
            end
            
            break
        end
        
        if norm(res) > maxval
            disp(['iterations = ', num2str(iter)]);
            error('Solution diverges');
        end
        
        if iter == iter_max
            disp([num2str(coordt(l))]);
            error('Does not converge');
        end
        
        iter = iter+1;
        T0 = T0 + dT; % update initial guess
        
    end
end
time=toc;

%%

for l=30:nt
    figure(1)
    hline = plot(coordx,T_final(:,l),'r','LineWidth',0.8);
    axis([xmin xmax 0 1100]);
    legend('numerical')
    title(sprintf('time = %d',coordt(l)))
    drawnow
    pause(0.2)
    M(M1) = getframe;
    M1 = M1+1;
end



