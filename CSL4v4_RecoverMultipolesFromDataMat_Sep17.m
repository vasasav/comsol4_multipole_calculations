% Written by Vassili Savinov on 29/06/12
% This function unscarbles the CMSL2_ComputeAllMultipoleData

% Modified on 17/07/12 to include the xiDip

% the first dimension is for different values for 
% different wavelengths

% Modified on 20/03/2013 to inculde electric/magnetic octupoles and
% toroidal dipole

% 15/06/2015
% changed title from RecoverMultipolesFromDataMat_CMSL3
% to CSL4v4_RecoverMultipolesFromDataMat
% to stick with the notation, no other changes

% 13/09/17 cheked by Vassili Savinov to bring it up to speed for lossy
% environment calculations

% Log30 p 116
function [MuPoles2D, MuPoles3D]=CSL4v4_RecoverMultipolesFromDataMat_Sep17(dataMat)
   % original encoding:
%    dataRow=...
%         [...
%         Id, Id_x, Id_y,  Im, Im_x, Im_y,  Idt, It, It_x,  It_y, IQ, IM,  Imr,...% 1-13
%         ...% 14-16
%         eDipVec_2D(1),      eDipVec_2D(2),      eDipVec_2D(3),...
%         ...% 17-19
%         mDipVec_2D(1),      mDipVec_2D(2),      mDipVec_2D(3),...
%         ...% 20-22
%         tDipVec_2D(1),      tDipVec_2D(2),      tDipVec_2D(3),...
%         ...% 23-25
%         rhoDipVec_2D(1),    rhoDipVec_2D(2),    rhoDipVec_2D(3),...
%         ...% 26-31
%         eQuadTens_2D(1,1),  eQuadTens_2D(2,2),  eQuadTens_2D(3,3), eQuadTens_2D(1,2), eQuadTens_2D(1,3), eQuadTens_2D(2,3),...
%         ...% 32-37
%         mQuadTens_2D(1,1),  mQuadTens_2D(2,2),  mQuadTens_2D(3,3), mQuadTens_2D(1,2), mQuadTens_2D(1,3), mQuadTens_2D(2,3),...
%         ...% 38-40
%         eDipVec_3D(1),      eDipVec_3D(2),      eDipVec_3D(3),...
%         ...% 41-43
%         mDipVec_3D(1),      mDipVec_3D(2),      mDipVec_3D(3),...
%         ...% 44-46
%         tDipVec_3D(1),      tDipVec_3D(2),      tDipVec_3D(3),...
%         ...% 47-49
%         rhoDipVec_3D(1),    rhoDipVec_3D(2),    rhoDipVec_3D(3),...
%         ...% 50-55
%         eQuadTens_3D(1,1),  eQuadTens_3D(2,2),  eQuadTens_3D(3,3), eQuadTens_3D(1,2), eQuadTens_3D(1,3), eQuadTens_3D(2,3),...
%         ...% 56-61
%         mQuadTens_3D(1,1),  mQuadTens_3D(2,2),  mQuadTens_3D(3,3), mQuadTens_3D(1,2), mQuadTens_3D(1,3), mQuadTens_3D(2,3),...
%         ... % 62-64 (here for historical reasons)
%         xiDip_2D(1),        xiDip_2D(2),        xiDip_2D(3),...
%         ... % 65-67
%         xiDip_3D(1),        xiDip_3D(2),        xiDip_3D(3),
% 
%         ... % 68-71
%         eOctXZZ_2D,     eOctYZZ_2D,     eOctXZZ_3D,     eOctYZZ_3D,...
%         ... % 72-75
%         tQuadXZ_2D,     tQuadYZ_2D,     tQuadXZ_3D,     tQuadYZ_3D,...
%         ... % 76-79
%         mOctXZZ_2D,     mOctYZZ_2D,     mOctXZZ_3D,     mOctYZZ_3D
%           ];
%% Multipoles from surface currents
    % electric dipole
    MuPoles2D.eDipX=dataMat(:,14);
    MuPoles2D.eDipY=dataMat(:,15);
    MuPoles2D.eDipZ=dataMat(:,16);
    
    % magnetic dipole
    MuPoles2D.mDipX=dataMat(:,17);
    MuPoles2D.mDipY=dataMat(:,18);
    MuPoles2D.mDipZ=dataMat(:,19);

    % toroidal dipole
    MuPoles2D.tDipX=dataMat(:,20);
    MuPoles2D.tDipY=dataMat(:,21);
    MuPoles2D.tDipZ=dataMat(:,22);
    
    % xi dipole
    MuPoles2D.xiDipX=dataMat(:,62);
    MuPoles2D.xiDipY=dataMat(:,63);
    MuPoles2D.xiDipZ=dataMat(:,64);
    
    % mean mag radius
    MuPoles2D.etaDipX=dataMat(:,23);
    MuPoles2D.etaDipY=dataMat(:,24);
    MuPoles2D.etaDipZ=dataMat(:,25);
    
    % electric quadrupole
    MuPoles2D.eQuadXX=dataMat(:,26);
    MuPoles2D.eQuadYY=dataMat(:,27);
    MuPoles2D.eQuadZZ=dataMat(:,28);
    MuPoles2D.eQuadXY=dataMat(:,29);
    MuPoles2D.eQuadXZ=dataMat(:,30);
    MuPoles2D.eQuadYZ=dataMat(:,31);
    
    % magnetic quadrupole
    MuPoles2D.mQuadXX=dataMat(:,32);
    MuPoles2D.mQuadYY=dataMat(:,33);
    MuPoles2D.mQuadZZ=dataMat(:,34);
    MuPoles2D.mQuadXY=dataMat(:,35);
    MuPoles2D.mQuadXZ=dataMat(:,36);
    MuPoles2D.mQuadYZ=dataMat(:,37);
    
    % electric octupole
    MuPoles2D.eOctXZZ=dataMat(:,68);
    MuPoles2D.eOctYZZ=dataMat(:,69);
    
    % toroidal quarupole
    MuPoles2D.tQuadXZ=dataMat(:,72);
    MuPoles2D.tQuadYZ=dataMat(:,73);
    
    % magnetic octupole
    MuPoles2D.mOctXZZ=dataMat(:,76);
    MuPoles2D.mOctYZZ=dataMat(:,77);
    
%% Multipoles from volume currents
    % electric dipole
    MuPoles3D.eDipX=dataMat(:,38);
    MuPoles3D.eDipY=dataMat(:,39);
    MuPoles3D.eDipZ=dataMat(:,40);
    
    % magnetic dipole
    MuPoles3D.mDipX=dataMat(:,41);
    MuPoles3D.mDipY=dataMat(:,42);
    MuPoles3D.mDipZ=dataMat(:,43);

    % toroidal dipole
    MuPoles3D.tDipX=dataMat(:,44);
    MuPoles3D.tDipY=dataMat(:,45);
    MuPoles3D.tDipZ=dataMat(:,46);
    
    % xi dipole
    MuPoles3D.xiDipX=dataMat(:,65);
    MuPoles3D.xiDipY=dataMat(:,66);
    MuPoles3D.xiDipZ=dataMat(:,67);
    
    % mean mag radius
    MuPoles3D.etaDipX=dataMat(:,47);
    MuPoles3D.etaDipY=dataMat(:,48);
    MuPoles3D.etaDipZ=dataMat(:,49);
    
    % electric quadrupole
    MuPoles3D.eQuadXX=dataMat(:,50);
    MuPoles3D.eQuadYY=dataMat(:,51);
    MuPoles3D.eQuadZZ=dataMat(:,52);
    MuPoles3D.eQuadXY=dataMat(:,53);
    MuPoles3D.eQuadXZ=dataMat(:,54);
    MuPoles3D.eQuadYZ=dataMat(:,55);
    
    % magnetic quadrupole
    MuPoles3D.mQuadXX=dataMat(:,56);
    MuPoles3D.mQuadYY=dataMat(:,57);
    MuPoles3D.mQuadZZ=dataMat(:,58);
    MuPoles3D.mQuadXY=dataMat(:,59);
    MuPoles3D.mQuadXZ=dataMat(:,60);
    MuPoles3D.mQuadYZ=dataMat(:,61);

    % electric octupole
    MuPoles3D.eOctXZZ=dataMat(:,70);
    MuPoles3D.eOctYZZ=dataMat(:,71);
    
    % toroidal quarupole
    MuPoles3D.tQuadXZ=dataMat(:,74);
    MuPoles3D.tQuadYZ=dataMat(:,75);
    
    % magnetic octupole
    MuPoles3D.mOctXZZ=dataMat(:,78);
    MuPoles3D.mOctYZZ=dataMat(:,79);
    
end