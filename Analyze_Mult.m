% Written by Vassili Savinov on 11/07/2019
% Analyze the multipole results for the mir-IR MM

clear;
close all;

%% load data
temp = dlmread('Bi2Te3_Mult_fixed_offset.csv');

% basics
lam_range = temp(:,1);
tran_amp_range = temp(:,2);% transmission amplitude
refl_amp_range = temp(:,3);% reflection amplitude

% incomplete multipoles
[~, incMuPoles]=CSL4v4_RecoverMultipolesFromDataMat_Sep17(temp(:,4:end));

% !!!!!!!!!!! In the origina simulation the wave is propagating along
% negative z, need to take it into account in calculating multipole
% reflection/transmission

%% fix multipoles
% the simulation in COMSOL is made in such a way that only half the
% structure, with positive y, is simulated. Therefore we need
% to complete the multipoles by flipping y-axis and adding the
% results for positive and negative y.
% this is a fairly non-standard thing to do, so I will do it on the case by
% case basis.

% I will only compute the multipoles that contribute to radiation here

% electric dipole, y-dipole vanishes, others double
MuPoles.eDipX       = 2 * incMuPoles.eDipX;
MuPoles.eDipY       = 0;

% magnetic dipole - easy way to remember is that when you flip
% the coordinate the Levi-Civita is multiplied by (-1) due to flipping 
% (usual tensor stuff) and then by another (-1) due to Jacobian change
% thus the net result of flipping y->-y is: +1 is for My -> My, 
% but -1 for Mx -> -Mx (and same for Z) since there only the Jacobian 
% sign change contributes
% this applies to all other multipoles
MuPoles.mDipX       = 0;
MuPoles.mDipY       = 2 * incMuPoles.mDipY;

% toroidal dipoles are polar, same as electric dipoles in terms of
% inversions
MuPoles.tDipX       = 2 * incMuPoles.tDipX;
MuPoles.tDipY       = 0;

% mean toroidal radius
MuPoles.xiDipX       = 2 * incMuPoles.xiDipX;
MuPoles.xiDipY       = 0;

% mean magnetic radius
MuPoles.etaDipX       = 0;
MuPoles.etaDipY       = 2 * incMuPoles.etaDipY;

% electric quadrupole. Trickier, polar in principle but two indices
% need to be considered separatelly
MuPoles.eQuadXZ       = 2 * incMuPoles.eQuadXZ;
MuPoles.eQuadYZ       = 0;

% magnetic quadrupole involves things like r_a (r x J)_b
% so when you fip y->-y:  
% _a -> _a if a!=y, and _a-> -_a if a=y
% and you get the opposite for _b 
MuPoles.mQuadXZ       = 0;  % ()()-> (+)(-) + (-)(+) 
MuPoles.mQuadYZ       = 2 * incMuPoles.mQuadYZ;  % ()()-> (+)(+) + (-)(-)

% toroidal quadrupole, like electric
MuPoles.tQuadXZ       = 2 * incMuPoles.tQuadXZ;
MuPoles.tQuadYZ       = 0;

% electric octupole
MuPoles.eOctXZZ       = 2 * incMuPoles.eOctXZZ;
MuPoles.eOctYZZ       = 0;

% magnetic octupole, two polar and one axial index
MuPoles.mOctXZZ       = 0;  % XZZ-> (+)(-)(+) 
MuPoles.mOctYZZ       = 2 * incMuPoles.mOctYZZ;  % ()()() -> (+)(-)(-) or (+)(+)(+)

% compute multipole fields. !!! NB at this point this should only be used
% to compre the relative power scattered by different multipoles
unit_area = 6.55e-6;
prop_dist = 0;

[MuPolField_X, MuPolField_Y,...
    ~, ~]=CSL4v4_ArbMed_GetAnalyticMuPolFields_Sep17(lam_range, unit_area, prop_dist, MuPoles);


%% visualize

figure(1);
ax=plotyy(lam_range/1e-6, abs(tran_amp_range).^2*100, ...
    lam_range/1e-6, abs(refl_amp_range).^2*100);
xlabel('Wavelength (um)');
ylabel(ax(1), 'Transmission (%)');
ylabel(ax(2), 'Reflection (%)');
legend('Transmission', 'Reflection');

lab_range={'eDip', 'mDip', 'tDip', 'eQuad', 'mQuad', 'tQuad', 'eOct', 'mOct'};
figure(2);
% subplot(121);
semilogy(lam_range/1e-6, ...
    abs([   MuPolField_X.eDip,  MuPolField_X.mDip,  MuPolField_X.tDip,  ...
        MuPolField_X.eQuad, MuPolField_X.mQuad, MuPolField_X.tQuad, ...
        MuPolField_X.eOct,  MuPolField_X.mOct]).^2 ...
    );
xlabel('Wavelength (um)');
ylabel('Power from arrays of multipoles. X-POLARIZATIOn');
legend(lab_range);
% %
% subplot(122);
% semilogy(lam_range/1e-6, ...
%     abs([   MuPolField_Y.eDip,  MuPolField_Y.mDip,  MuPolField_Y.tDip,  ...
%         MuPolField_Y.eQuad, MuPolField_Y.mQuad, MuPolField_Y.tQuad, ...
%         MuPolField_Y.eOct,  MuPolField_Y.mOct]).^2 ...
%     );
% xlabel('Wavelength (um)');
% ylabel('Power from arrays of multipoles. Y-POLARIZATIOn');
% legend(lab_range);

x_minmax = [4e-6, 12e-6];
good_x_inds = logical( (lam_range>=x_minmax(1)) & (lam_range<=x_minmax(2)) );

y_mat=  abs([   MuPolField_X.eDip,  MuPolField_X.mDip,  MuPolField_X.tDip,  ...
        MuPolField_X.eQuad, MuPolField_X.mQuad, MuPolField_X.tQuad, ...
        MuPolField_X.eOct,  MuPolField_X.mOct]).^2;

y_max = max(max(y_mat(good_x_inds, :)));
y_minmax = [y_max*1e-6, y_max];

figure(11);
semilogy(lam_range/1e-6, y_mat/y_max);
xlabel('Wavelength (um)');
ylabel('Power from arrays of multipoles. X-POLARIZATIOn');
legend(lab_range);
xlim(x_minmax/1e-6);
ylim(y_minmax/y_max);

% check phase between electric and toroidal dipoles
figure(12);
plot(lam_range/1e-6, angle(MuPolField_X.eDip./MuPolField_X.tDip)/pi*180);
ylabel('phase difference between electric and toroidal field emission (deg)');
xlabel('Wavelength (um)');
xlim(x_minmax/1e-6);
