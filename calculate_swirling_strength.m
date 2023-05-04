clear,clc

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));

ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));

load('velocity_data.mat');
nx = 112;
ny = 112;
nz = 128;
dkx = 1;
dky = 2;
Lx = 2*pi;
Ly = 1*pi;
Retau = 180;
xgrid = linspace(0, Lx, nx+1);
ygrid = linspace(0, Ly, ny+1);
xp = xgrid(1:nx)/2 + xgrid(2:nx+1)/2;
yp = ygrid(1:ny)/2 + ygrid(2:ny+1)/2; 

[Diff,zc] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

kx = [0:dkx:(50*dkx),-(50*dkx):dkx:-dkx];
ky = [0:dky:(50*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

[velocity_tensor] = real(get_velocity_tensor(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky));
velocity_tensor = permute(velocity_tensor,[4,1,2,3]); 
[swirling_strength] = get_swirling_strength(velocity_tensor,nx,ny,nz);

swirling_strength = permute(swirling_strength,[2,3,1]); %[y,x,z]
u = permute(u,[2,3,1]);
v = permute(v,[2,3,1]);
w = permute(w,[2,3,1]);

x_grid_array = xp;
y_grid_array = yp;
z_grid_array = zc(2:end-1);


%%
max(swirling_strength,[],'all')
figure 
s = isosurface(x_grid_array,y_grid_array,z_grid_array,swirling_strength,0.33*max(swirling_strength,[],'all')); hold on
p = patch(s);
isonormals(x_grid_array,y_grid_array,z_grid_array,swirling_strength,p)
view(3);
set(p,'FaceColor','y');  
set(p,'EdgeColor','none');
lightangle(gca,80,20)
lighting gouraud
axis([x_grid_array(1) x_grid_array(end) y_grid_array(1) y_grid_array(end) z_grid_array(end) z_grid_array(1)])
daspect([1 1 1])
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca,'Fontsize',16)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';