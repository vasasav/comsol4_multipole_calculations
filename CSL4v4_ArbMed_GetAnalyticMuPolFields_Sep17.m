% Written by Vassili Savinov on 09/09/12
% Compute the back- and forward-scattered fields
% of the MM from  the multipoles
% using the analytic relations

% NB! back-scattered fields relate to 
% reflection whilst the forward scattered fields relate to transmission

% Assuming that incident radiation
% is the plane wave propagating along the positive z-axis

% unitArea - is the area of the unit cell of the metamaterial
% lambda is the wavelength 

% observer-dist: how far is the observer - needed to 
% get the correct phase

% Details in LogBook 10 see pages 66, 67, 69, 70 and related

% Updated on 20/03/2013 to include electric/magnetic octupoles
% and the toroidal quadrupole see LogBook 12 p 38, 43, 51, 92

% Name changed from GetAnalyticMuPolFields_CMSL3 on 16/06/15
% No Other changes

% 20/04/17 extended to handle calculations for arbitary medium surrounding
% the metamaterial.

% 13/09/17 cheked by Vassili Savinov to bring it up to speed for lossy
% environment calculations

% Log30 p 116

function [MuPolForwF_Xpol, MuPolForwF_Ypol,...
    MuPolBackF_Xpol, MuPolBackF_Ypol]=CSL4v4_ArbMed_GetAnalyticMuPolFields_Sep17(lambda, unitArea, observerDist, MuPoles, eps1)

    if ~exist('eps1', 'var')
        eps1=1;% back to vacuum
    end

    mu0=4*pi*1e-7;
    lightC=2.99792458e8/sqrt(eps1); % Sep17 changed to exact value, compensate for envoronment
    wavn=2*pi*sqrt(eps1)./lambda;% wavenumber
    % the unit-setting constant and the correct phase
    % due to propgation towards the obesrcer
    preConst=mu0*lightC^2/(2*unitArea).*exp(-1i*wavn.*observerDist);
    
    %% X-polarization
    % firstly compute the fields for 
    % forward scattering LogBook 10 p 66
    RdotZ=+1;% forward scattering
    % electric dipole
    MuPolForwF_Xpol.eDip=preConst.*(-1i*wavn).*(MuPoles.eDipX);
    % magnetic dipole
    MuPolForwF_Xpol.mDip=preConst.*(-1i*wavn).*...
        (MuPoles.mDipY-(wavn.^2/10).*MuPoles.etaDipY)*RdotZ;
    % toroidal dipole
    MuPolForwF_Xpol.tDip=preConst.*(-wavn.^2).*...
        (MuPoles.tDipX+(wavn.^2/10).*MuPoles.xiDipX);
    % electric quadrupole
    MuPolForwF_Xpol.eQuad=preConst.*(wavn.^2).*(MuPoles.eQuadXZ)*RdotZ;

    % magnetic quadrupole
    MuPolForwF_Xpol.mQuad=preConst.*(wavn.^2/2).*(MuPoles.mQuadYZ);
    
    % electric octupole
    MuPolForwF_Xpol.eOct=preConst.*(1i*wavn.^3).*MuPoles.eOctXZZ;
    
    % toroidal quadrupole
    MuPolForwF_Xpol.tQuad=preConst.*(-1i*wavn.^3/3).*(MuPoles.tQuadXZ)*RdotZ;

    % magnetic octupole
    MuPolForwF_Xpol.mOct=preConst.*(-1i*wavn.^3/180).*(-MuPoles.mOctYZZ)*RdotZ;
    
    MuPolForwF_Xpol.all_mult=   MuPolForwF_Xpol.eDip+...
                                MuPolForwF_Xpol.mDip+...
                                MuPolForwF_Xpol.tDip+...
                                MuPolForwF_Xpol.eQuad+...
                                MuPolForwF_Xpol.mQuad+...
                                MuPolForwF_Xpol.eOct+...
                                MuPolForwF_Xpol.tQuad+...
                                MuPolForwF_Xpol.mOct;
    
    % for reverse scattering simply reverse the RdotZ
    % I will not recompute it, but will simply reverse the signs
    RdotZ=-1;% backward scattering
    % electric dipole
    MuPolBackF_Xpol.eDip=MuPolForwF_Xpol.eDip;
    % magnetic dipole
    MuPolBackF_Xpol.mDip=MuPolForwF_Xpol.mDip*RdotZ;
    % toroidal dipole
    MuPolBackF_Xpol.tDip=MuPolForwF_Xpol.tDip;
    % electric quadrupole
    MuPolBackF_Xpol.eQuad=MuPolForwF_Xpol.eQuad*RdotZ;
    % magnetic quadrupole
    MuPolBackF_Xpol.mQuad=MuPolForwF_Xpol.mQuad;
    
    % electric octupole
    MuPolBackF_Xpol.eOct=preConst.*(1i*wavn.^3).*MuPoles.eOctXZZ;
    
    % toroidal quadrupole
    MuPolBackF_Xpol.tQuad=preConst.*(-1i*wavn.^3/3).*(MuPoles.tQuadXZ)*RdotZ;

    % magnetic octupole
    MuPolBackF_Xpol.mOct=preConst.*(-1i*wavn.^3/180).*(-MuPoles.mOctYZZ)*RdotZ;

    MuPolBackF_Xpol.all_mult=   MuPolBackF_Xpol.eDip+...
                                MuPolBackF_Xpol.mDip+...
                                MuPolBackF_Xpol.tDip+...
                                MuPolBackF_Xpol.eQuad+...
                                MuPolBackF_Xpol.mQuad+...
                                MuPolBackF_Xpol.eOct+...
                                MuPolBackF_Xpol.tQuad+...
                                MuPolBackF_Xpol.mOct;
    
    %% Y-polarization
    % firstly compute the fields for 
    % forward scattering LogBook 10 p 69
    RdotZ=+1;% forward scattering
    % electric dipole
    MuPolForwF_Ypol.eDip=preConst.*(-1i*wavn).*(MuPoles.eDipY);
    % magnetic dipole
    MuPolForwF_Ypol.mDip=preConst.*(+1i*wavn).*...
        (MuPoles.mDipX-(wavn.^2/10).*MuPoles.etaDipX)*RdotZ;
    % toroidal dipole
    MuPolForwF_Ypol.tDip=preConst.*(-wavn.^2).*...
        (MuPoles.tDipY+(wavn.^2/10).*MuPoles.xiDipY);
    % electric quadrupole
    MuPolForwF_Ypol.eQuad=preConst.*(wavn.^2).*(MuPoles.eQuadYZ)*RdotZ;
    % magnetic quadrupole
    MuPolForwF_Ypol.mQuad=preConst.*(-wavn.^2/2).*(MuPoles.mQuadXZ);
    
    % electric octupole
    MuPolForwF_Ypol.eOct=preConst.*(1i*wavn.^3).*MuPoles.eOctYZZ;

    % toroidal quadrupole
    MuPolForwF_Ypol.tQuad=preConst.*(-1i*wavn.^3/3).*(MuPoles.tQuadYZ)*RdotZ;

    % magnetic octupole
    MuPolForwF_Ypol.mOct=preConst.*(-1i*wavn.^3/180).*(+MuPoles.mOctXZZ)*RdotZ;

    MuPolForwF_Ypol.all_mult=   MuPolForwF_Ypol.eDip+...
                                MuPolForwF_Ypol.mDip+...
                                MuPolForwF_Ypol.tDip+...
                                MuPolForwF_Ypol.eQuad+...
                                MuPolForwF_Ypol.mQuad+...
                                MuPolForwF_Ypol.eOct+...
                                MuPolForwF_Ypol.tQuad+...
                                MuPolForwF_Ypol.mOct;
    
    % for reverse scattering simply reverse the RdotZ
    % I will not recompute it, but will simply reverse the signs
    RdotZ=-1;% backward scattering
    % electric dipole
    MuPolBackF_Ypol.eDip=MuPolForwF_Ypol.eDip;
    % magnetic dipole
    MuPolBackF_Ypol.mDip=MuPolForwF_Ypol.mDip*RdotZ;
    % toroidal dipole
    MuPolBackF_Ypol.tDip=MuPolForwF_Ypol.tDip;
    % electric quadrupole
    MuPolBackF_Ypol.eQuad=MuPolForwF_Ypol.eQuad*RdotZ;
    % magnetic quadrupole
    MuPolBackF_Ypol.mQuad=MuPolForwF_Ypol.mQuad;

    % electric octupole
    MuPolBackF_Ypol.eOct=preConst.*(1i*wavn.^3).*MuPoles.eOctYZZ;

    % toroidal quadrupole
    MuPolBackF_Ypol.tQuad=preConst.*(-1i*wavn.^3/3).*(MuPoles.tQuadYZ)*RdotZ;

    % magnetic octupole
    MuPolBackF_Ypol.mOct=preConst.*(-1i*wavn.^3/180).*(+MuPoles.mOctXZZ)*RdotZ;
    
    MuPolBackF_Ypol.all_mult=   MuPolBackF_Ypol.eDip+...
                                MuPolBackF_Ypol.mDip+...
                                MuPolBackF_Ypol.tDip+...
                                MuPolBackF_Ypol.eQuad+...
                                MuPolBackF_Ypol.mQuad+...
                                MuPolBackF_Ypol.eOct+...
                                MuPolBackF_Ypol.tQuad+...
                                MuPolBackF_Ypol.mOct;
end