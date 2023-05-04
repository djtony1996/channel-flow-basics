clear,clc
fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));

load('velocity_data.mat');
nx = 112;
ny = 112;
nz = 128;
dkx = 1;
dky = 2;
Retau = 180;

[Diff,~] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);
[~,WEIGHT] = clenCurt(nz);
WEIGHT = WEIGHT(2:nz);

k_edge = 50;
kx = [0:dkx:(k_edge*dkx),-(k_edge*dkx):dkx:-dkx];
ky = [0:dky:(k_edge*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

kx_0posi = [0:dkx:20*dkx];
ky_0posi = [0:dky:20*dky];
kx_posi  = kx_0posi(2:end);
ky_posi  = ky_0posi(2:end);

u_F  = fft_xy(u);
v_F  = fft_xy(v);
w_F  = fft_xy(w);

[duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
[dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
[dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m); 

M_4d     = zeros(length(ky_0posi),length(kx_0posi),length(ky_0posi),length(kx_0posi));

for kx_wave = 1: length(kx_0posi)
    for ky_wave = 1: length(ky_0posi)
        k_x = kx_0posi(kx_wave);
        k_y = ky_0posi(ky_wave);

        if k_x == 0 && k_y == 0
           continue
        end

        if k_x == 0
            M_single = (-2) .* real(get_N_transfer_0b(k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
        elseif k_y == 0
            M_single = (-2) .* real(get_N_transfer_a0(k_x,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
        else 
            M_single = (-2) .* real(get_N_transfer_ab(k_x,k_y,u_F,v_F,w_F,duF_dx,duF_dy,duF_dz,dvF_dx,dvF_dy,dvF_dz,dwF_dx,dwF_dy,dwF_dz,kx_posi,ky_posi,nx,ny,dkx,dky,WEIGHT));
        end    
        M_4d(:,:,ky_wave,kx_wave) = M_single;
    end
end
