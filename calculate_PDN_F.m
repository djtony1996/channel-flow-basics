clear,clc

load('velocity_data.mat');
nx = 112;
ny = 112;
nz = 128;
dkx = 1;
dky = 2;
Retau = 180;

[Diff,~] = cheb(nz);
Diff = Diff(2:end-1,2:end-1);

k_edge = 50;
kx = [0:dkx:(k_edge*dkx),-(k_edge*dkx):dkx:-dkx];
ky = [0:dky:(k_edge*dky)];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

% calculate the energy in Fourier space
[ProF, DissF, NonTF, ~, ~, ~] = get_three_energy(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky,Retau,dU_mean_dz,Diff);
NonTF = real(NonTF); 

temp = ProF(1:k_edge,1:k_edge);
ProF(2:k_edge,2:k_edge)  = temp(2:k_edge,2:k_edge) + fliplr(ProF(2:k_edge,(end-k_edge+2):end));
temp = DissF(1:k_edge,1:k_edge);
DissF(2:k_edge,2:k_edge) = temp(2:k_edge,2:k_edge) + fliplr(DissF(2:k_edge,(end-k_edge+2):end));
temp = NonTF(1:k_edge,1:k_edge);
NonTF(2:k_edge,2:k_edge) = temp(2:k_edge,2:k_edge) + fliplr(NonTF(2:k_edge,(end-k_edge+2):end));

