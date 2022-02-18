% Written by Vassili Savinov on 01/03/12
% This function parses the strings into
% form which COMSOL will understand

% Upgraded to parse for COMSOL 4.4 on 11/06/15
% by Vassili Savinov
% for now will only implement the volume current densities
% since the suraface current densities are a known pest

% convention
% e0 = epsilon0_const
% m0 = mu0_const
% light= c_const

% will only provide for isotropic materials for now
% erX = emw.epsilonrxx
% erY = emw.epsilonryy
% erZ = emw.epsilonrzz
% omg = emw.omega

% JVx, JVy, JVz = i*omg*ep0*(eprX-1)*emw.Ex, 
% x,y,z = [untouched]
% r2 = (x^2+y^2+z^2)
% R.J = J.R = (x*JVx+y*JVy+z*JVz)
% RxJ_x = -JxR_x = (y*JVz-z*JVy), ...

% 20/04/17 extended to handle calculations for arbitary medium surrounding
% the metamaterial. In essence this is to allow free, nonmetamaterial space
% to be not vacuum. To do this I simulatneously modify the definition of
% the current density and the speed of light
% in essence, COMSOL solves curl(B)=(i\omega\epsilon/c^2)E
% I split it into curl(B)=\mu_0*(i\omega\epsilon_0*(espilon-eps1)E)+i\omega/(c/sqrt(eps1))^2 E
% therefore the current density is now J=i\omega\epsilon_0*(espilon-eps1)E
% and speed of light v=c/sqrt(eps1) where eps1 is the background dielectic
% constant
% see Log 29.68

% this should reduce the likelehood of typing the wrong formula
% eps1 is the background dielectric constant

% 13/09/17 cheked by Vassili Savinov to bring it up to speed for lossy
% environment calculations

% Log30 p 116

function retStr=CSL4v4_ArbMed_GetMultipoleIntegrand_Sep17(parStr, eps1)
    if ~exist('eps1', 'var')
        eps1=1;% back to vacuum
    end

    iStart=1000; % starting position of the found token
    
    % prepare the string 
    retStr=parStr;
    
    %% prepare the dictionary
    findColl={};% find collection
    replColl={};% replace collection
    
    iEntry=0;
    % populate the dictionary
    % note that it is self-referencible 
    % i.e. you can use definitions of dictionary 
    % within other definitions
    iEntry=iEntry+1;
    findColl{iEntry}='e0'; % what is taken out
    replColl{iEntry}='epsilon0_const'; % what is put in instead (for COMSOl)
    
    iEntry=iEntry+1;
    findColl{iEntry}='m0';
    replColl{iEntry}='mu0_const';
    
    iEntry=iEntry+1;
    findColl{iEntry}='light';
%     replColl{iEntry}='c_const';
    replColl{iEntry}=sprintf('(c_const/(%.6f+i*(%.6f)))', real(sqrt(eps1)), imag(sqrt(eps1)));
    
    iEntry=iEntry+1;
    findColl{iEntry}='erX';
    replColl{iEntry}='emw.epsilonrxx';
    
    iEntry=iEntry+1;
    findColl{iEntry}='erY';
    replColl{iEntry}='emw.epsilonryy';
    
    iEntry=iEntry+1;
    findColl{iEntry}='erZ';
    replColl{iEntry}='emw.epsilonrzz';
    
    iEntry=iEntry+1;
    findColl{iEntry}='omg';
    replColl{iEntry}='emw.omega';
    
    % this definition of volume current density
    % only uses the comlex dielectric constant 
    % of the material
    % it should be possible to extend it to 
    % a more general case of ((i*omg*e0*(epr-1)+sigma_rfw)*Ex
    % but I have no use for it now, so I'll update the function when the
    % need will arrise (and will simultaneously test this upgraded functionality)
    iEntry=iEntry+1;
    findColl{iEntry}='JVx';
    replColl{iEntry}=sprintf('(i*omg*e0*(erX-(%.6f+i*(%.6f)))*emw.Ex)', real(eps1), imag(eps1));
    
    iEntry=iEntry+1;
    findColl{iEntry}='JVy';
    replColl{iEntry}=sprintf('(i*omg*e0*(erY-(%.6f+i*(%.6f)))*emw.Ey)', real(eps1), imag(eps1));
    
    iEntry=iEntry+1;
    findColl{iEntry}='JVz';
    replColl{iEntry}=sprintf('(i*omg*e0*(erZ-(%.6f+i*(%.6f)))*emw.Ez)', real(eps1), imag(eps1));
    
    iEntry=iEntry+1;
    findColl{iEntry}='R2';
    replColl{iEntry}='(x^2+y^2+z^2)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='R.JV';
    replColl{iEntry}='(x*JVx+y*JVy+z*JVz)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='JV.R';
    replColl{iEntry}='R.JV';
    
    iEntry=iEntry+1;
    findColl{iEntry}='RxJV_x';
    replColl{iEntry}='(y*JVz-z*JVy)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='RxJV_y';
    replColl{iEntry}='(z*JVx-x*JVz)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='RxJV_z';
    replColl{iEntry}='(x*JVy-y*JVx)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='JVxR_x';
    replColl{iEntry}='(-RxJV_x)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='JVxR_y';
    replColl{iEntry}='(-RxJV_y)';
    
    iEntry=iEntry+1;
    findColl{iEntry}='JVxR_z';
    replColl{iEntry}='(-RxJV_z)';
    
    %% evaluate
     
    while ~isempty(iStart)
        
        % check the string for key phrases. go through the collection
        for iEntry=1:length(findColl)
            iStart=strfind(retStr, findColl{iEntry});
            
            % found a match
            if ~isempty(iStart)
                iStart=iStart(1);% concentrate on the first occurance
                % only, other ones will be analyzed next time 
                
                % replace this occurance from the collection
                entLen=length(findColl{iEntry});
                retStr=sprintf('%s%s%s', retStr(1:(iStart-1)), replColl{iEntry},...
                retStr((iStart+entLen):end));
                
                break;
            end
        end
    end
end