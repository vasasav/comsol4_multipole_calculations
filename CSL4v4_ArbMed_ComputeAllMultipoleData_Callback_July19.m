% Written by Vassili Savinov on 03/03/12
% compute all multipoles and all intensities using COMSOL 
% according to formulae derrived by me
% integrates both the 2D and the 3d multipoles 
% uses GetMultipoleIntegrand to produce long integrand 
% expression for COMSOL from easy-to-check expressions I enter

% upgraded on 19/07/12 to include xi dipole
% rho is replaced by eta to comply with new convention

% upgrades on 20/03/13 to include electric octupole, magnetic octupole
% and toroidal quadrupole

% updated on 20/04/15 eQuadZZ error has been corrected 
% (used to be asigned into YY entry)

% Note that one stray version of the CMSL3 code also 
% had non-assignment of mOctXZZ_3D, however this bug did not 
% make it to PRB paper (4 pairs of double SRRs) or to 
% planar superconducting MM

% Upgraded on 11/06/2015 to run on COMSOL 4.4
% will only support volume currents for now

% 20/04/17 extended to handle calculations for arbitary medium surrounding
% the metamaterial.

% updated on 04/05/17 
% now the data is assigned correctly into [eQuadTens_3D(2,2) entry 
% Thanks Wei-Yi!

% 13/09/17 cheked by Vassili Savinov to bring it up to speed for lossy
% environment calculations

% Log30 p 116

% June' 19 corrected the equations for calculating electric otcupole

% July 10/07/2019 I am tired of adapting to the 
% perculiarities of each comsol model. Will now do integration via
% call-back [val, unit]=my_int_callback('string')

function dataRow=CSL4v4_ArbMed_ComputeAllMultipoleData_Callback_July19(integration_callback, fq, eps1)

    % all the code for 2d multipoles is kept for the sake of
    % backwards-compat

    if ~exist('eps1', 'var')
        eps1=1;% back to vacuum
    end

    lightC=2.99792458e8/sqrt(eps1);%Sep17
    om=2*pi*fq;
    
    eDipVec_2D=zeros(1,3);      eDipVec_3D=zeros(1,3);
    eDipUnit='C*m';
    
    mDipVec_2D=zeros(1,3);      mDipVec_3D=zeros(1,3);
    mDipUnit='C*m';
    
    tDipVec_2D=zeros(1,3);      tDipVec_3D=zeros(1,3);
    tDipUnit='m^2*s*A';
    
    xiDipVec_2D=zeros(1,3);     xiDipVec_3D=zeros(1,3);
    xiDipUnit='m^4*s*A';
    
    etaDipVec_2D=zeros(1,3);    etaDipVec_3D=zeros(1,3);
    etaDipUnit='m^3*s*A';
    
    eQuadTens_2D=zeros(3,3);    eQuadTens_3D=zeros(3,3);
    eQuadTensUnit='m^2*s*A';
    
    mQuadTens_2D=zeros(3,3);    mQuadTens_3D=zeros(3,3);
    mQuadTensUnit='m^2*s*A';

    eOctXZZ_2D=0;
    eOctYZZ_2D=0;
    eOctXZZ_3D=0;
    eOctYZZ_3D=0;
    eOctTensUnit='m^3*s*A';
    
    tQuadXZ_2D=0;
    tQuadYZ_2D=0;
    tQuadXZ_3D=0;
    tQuadYZ_3D=0;
    tQuadTensUnit='m^3*s*A';
    
    mOctXZZ_2D=0;
    mOctYZZ_2D=0;
    mOctXZZ_3D=0;
    mOctYZZ_3D=0;
    mOctTensUnit='m^3*s*A';
    
    bad_units = false;
    
    % unlike in COMSOL 3 one cannot choose units of integration in
    % COMSOL 4, we can just check them

    % electric dipole
    [eDipVec_3D(1), unit]=integration_callback( CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('JVx/(i*omg)',eps1) );                    bad_units = bad_units || ~strcmp(unit,eDipUnit);
    [eDipVec_3D(2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('JVy/(i*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eDipUnit);
    [eDipVec_3D(3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('JVz/(i*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eDipUnit);

    % magnetic dipole
    [mDipVec_3D(1), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('RxJV_x/(2*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mDipUnit);
    [mDipVec_3D(2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('RxJV_y/(2*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mDipUnit);
    [mDipVec_3D(3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('RxJV_z/(2*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mDipUnit);

    % toroidal dipole
    [tDipVec_3D(1), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(R.JV*x-2*R2*JVx)/(10*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,tDipUnit);
    [tDipVec_3D(2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(R.JV*y-2*R2*JVy)/(10*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,tDipUnit); 
    [tDipVec_3D(3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(R.JV*z-2*R2*JVz)/(10*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,tDipUnit);

    % mean toroidal radius
    [xiDipVec_3D(1), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('R2*(3*R2*JVx-2*x*R.JV)/(28*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,xiDipUnit);
    [xiDipVec_3D(2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('R2*(3*R2*JVy-2*y*R.JV)/(28*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,xiDipUnit);
    [xiDipVec_3D(3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('R2*(3*R2*JVz-2*z*R.JV)/(28*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,xiDipUnit);

    % mean magnetic radius
    [etaDipVec_3D(1), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(RxJV_x*R2)/(2*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,etaDipUnit);
    [etaDipVec_3D(2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(RxJV_y*R2)/(2*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,etaDipUnit);
    [etaDipVec_3D(3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(RxJV_z*R2)/(2*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,etaDipUnit);
    % electric quadrupole
    [eQuadTens_3D(1,1), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(x*JVx-R.JV/3)/(i*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eQuadTensUnit);
    [eQuadTens_3D(2,2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(y*JVy-R.JV/3)/(i*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eQuadTensUnit);
    [eQuadTens_3D(3,3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(z*JVz-R.JV/3)/(i*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eQuadTensUnit);
    %
    [eQuadTens_3D(1,2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(x*JVy+y*JVx)/(i*2*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eQuadTensUnit);
    [eQuadTens_3D(1,3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(x*JVz+z*JVx)/(i*2*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eQuadTensUnit);
    [eQuadTens_3D(2,3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(y*JVz+z*JVy)/(i*2*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eQuadTensUnit);
    eQuadTens_3D(2,1)=eQuadTens_3D(1,2);
    eQuadTens_3D(3,1)=eQuadTens_3D(1,3);
    eQuadTens_3D(3,2)=eQuadTens_3D(2,3);

    % magnetic quadrupole
    [mQuadTens_3D(1,1), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(2*RxJV_x*x)/(3*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mQuadTensUnit);
    [mQuadTens_3D(2,2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(2*RxJV_y*y)/(3*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mQuadTensUnit);
    [mQuadTens_3D(3,3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(2*RxJV_z*z)/(3*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mQuadTensUnit);
    %
    [mQuadTens_3D(1,2), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(RxJV_x*y+x*RxJV_y)/(3*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mQuadTensUnit);
    [mQuadTens_3D(1,3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(RxJV_x*z+x*RxJV_z)/(3*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mQuadTensUnit);
    [mQuadTens_3D(2,3), unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('(RxJV_y*z+y*RxJV_z)/(3*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,mQuadTensUnit);
    %
    mQuadTens_3D(2,1)=mQuadTens_3D(1,2);
    mQuadTens_3D(3,1)=mQuadTens_3D(1,3);
    mQuadTens_3D(3,2)=mQuadTens_3D(2,3);

    % electric octupole, only the relevant compnents are found
    [eOctXZZ_3D, unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('( (z^2-R2/5)*JVx + (2*x*z)*JVz - (2/5)*(x+2*z)*R.JV )/(i*6*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eOctTensUnit);
    [eOctYZZ_3D, unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('( (z^2-R2/5)*JVy + (2*y*z)*JVz - (2/5)*(y+2*z)*R.JV )/(i*6*omg)',eps1) );                   bad_units = bad_units || ~strcmp(unit,eOctTensUnit);
    % toroidal quadrupole
    [tQuadXZ_3D, unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('( (4*x^2*z-5*z*R2)*JVx + (4*x*y*z)*JVy + (4*x*z^2-5*x*R2)*JVz )/(28*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,tQuadTensUnit);
    [tQuadYZ_3D, unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('( (4*x*y*z)*JVx + (4*y^2*z-5*z*R2)*JVy + (4*y*z^2-5*y*R2)*JVz )/(28*light)',eps1) );                   bad_units = bad_units || ~strcmp(unit,tQuadTensUnit);

    % magnetic octupole
    [mOctXZZ_3D, unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('( (-10*x*y*z)*JVx + (11*x^2*z+y^2*z-4*z^3)*JVy + (4*y*z^2-x^2*y-y^3)*JVz )*(3/(2*light))',eps1) );                   bad_units = bad_units || ~strcmp(unit,mOctTensUnit);
    [mOctYZZ_3D, unit]=integration_callback(  CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17('( (4*z^3-x^2*z-11*y^2*z)*JVx + (10*x*y*z)*JVy + (x^3+x*y^2-4*x*z^2)*JVz )*(3/(2*light))',eps1) );                   bad_units = bad_units || ~strcmp(unit,mOctTensUnit);

    
    if bad_units
        error('CSL4v4_ComputeAllMultipoleData: Wrong units!');
    end

    %%%
    
    % total multipoles
    eDipVec=eDipVec_2D+eDipVec_3D;
    mDipVec=mDipVec_2D+mDipVec_3D;
    tDipVec=tDipVec_2D+tDipVec_3D;
    etaDipVec=etaDipVec_2D+etaDipVec_3D;
    eQuadTens=eQuadTens_2D+eQuadTens_3D;
    mQuadTens=mQuadTens_2D+mQuadTens_3D;
    
    % now compute each intensity term
    Id=(2*om^4)/(3*lightC^3)*sum(abs(eDipVec).^2);
    Im=(2*om^4)/(3*lightC^3)*sum(abs(mDipVec).^2);
    Idt=(4*om^5)/(3*lightC^4)*imag(sum(conj(eDipVec).*tDipVec));
    It=(2*om^6)/(3*lightC^5)*sum(abs(tDipVec).^2);
    IQ=(om^6)/(5*lightC^5)*sum(sum(abs(eQuadTens).^2));    
    IM=(om^6)/(20*lightC^5)*sum(sum(abs(mQuadTens).^2));
    Imr=-(2*om^6)/(15*lightC^5)*real(sum(conj(mDipVec).*etaDipVec));
    
    % extra terms resolved into components
    Id_x=(2*om^4)/(3*lightC^3)*abs(eDipVec(1))^2;
    Id_y=(2*om^4)/(3*lightC^3)*abs(eDipVec(2))^2;
    
    Im_x=(2*om^4)/(3*lightC^3)*abs(mDipVec(1))^2;
    Im_y=(2*om^4)/(3*lightC^3)*abs(mDipVec(2))^2;
    
    It_x=(2*om^6)/(3*lightC^5)*abs(tDipVec(1))^2;
    It_y=(2*om^6)/(3*lightC^5)*abs(tDipVec(2))^2;
    
    % row containing all the data
    % convinient for saving in file 
    dataRow=...
        [...
        Id, Id_x, Id_y,  Im, Im_x, Im_y,  Idt, It, It_x,  It_y, IQ, IM,  Imr,...% 1-13
        ...% 14-19
        eDipVec_2D(1),      eDipVec_2D(2),      eDipVec_2D(3),...
        mDipVec_2D(1),      mDipVec_2D(2),      mDipVec_2D(3),...
        ...% 20-37
        tDipVec_2D(1),      tDipVec_2D(2),      tDipVec_2D(3),...
        etaDipVec_2D(1),    etaDipVec_2D(2),    etaDipVec_2D(3),...
        eQuadTens_2D(1,1),  eQuadTens_2D(2,2),  eQuadTens_2D(3,3), eQuadTens_2D(1,2), eQuadTens_2D(1,3), eQuadTens_2D(2,3),...
        mQuadTens_2D(1,1),  mQuadTens_2D(2,2),  mQuadTens_2D(3,3), mQuadTens_2D(1,2), mQuadTens_2D(1,3), mQuadTens_2D(2,3),...
        ...% 38-61
        eDipVec_3D(1),      eDipVec_3D(2),      eDipVec_3D(3),...
        mDipVec_3D(1),      mDipVec_3D(2),      mDipVec_3D(3),...
        tDipVec_3D(1),      tDipVec_3D(2),      tDipVec_3D(3),...
        etaDipVec_3D(1),    etaDipVec_3D(2),    etaDipVec_3D(3),...
        eQuadTens_3D(1,1),  eQuadTens_3D(2,2),  eQuadTens_3D(3,3), eQuadTens_3D(1,2), eQuadTens_3D(1,3), eQuadTens_3D(2,3),...
        mQuadTens_3D(1,1),  mQuadTens_3D(2,2),  mQuadTens_3D(3,3), mQuadTens_3D(1,2), mQuadTens_3D(1,3), mQuadTens_3D(2,3),...
        ... % 62-64 (here for historical reasons)
        xiDipVec_2D(1),     xiDipVec_2D(2),     xiDipVec_2D(3),...
        ... % 65-67
        xiDipVec_3D(1),     xiDipVec_3D(2),     xiDipVec_3D(3)...
        ... % 68-71
        eOctXZZ_2D,     eOctYZZ_2D,     eOctXZZ_3D,     eOctYZZ_3D,...
        ... % 72-75
        tQuadXZ_2D,     tQuadYZ_2D,     tQuadXZ_3D,     tQuadYZ_3D,...
        ... % 76-79
        mOctXZZ_2D,     mOctYZZ_2D,     mOctXZZ_3D,     mOctYZZ_3D
        ];
end