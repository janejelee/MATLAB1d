
n=50;
phi = 0.25;

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

valrhof = interp1(rhof_raw(:,1),rhof_raw(:,2),linspace(0,1000,n));
valcf = interp1(cf_raw(:,1),cf_raw(:,2),linspace(0,1000,n));
valcr = interp1(cr_raw(:,1),cr_raw(:,2),linspace(0,1000,n));
valkappaf = interp1(kappaf_raw(:,1),kappaf_raw(:,2),linspace(0,1000,n));
valkappar = interp1(kappar_raw(:,1),kappar_raw(:,2),linspace(0,1000,n));

kappa_val = valkappaf.^phi .* valkappar.^(1-phi);

coord = linspace(0,10000,n);

%%

coord= coord';
valrhof = valrhof';
kappa_val = kappa_val';
valcf = valcf';
valcr = valcr';


f_kappa = fit(coord,kappa_val,'exp2');
f_rhof = fit(coord,valrhof,'exp2');
f_cf = fit(coord,valcf,'exp2');
f_cr = fit(coord,valcr,'exp2');
plot(f,coord,kappa_val)



