function [  M_1_V, Tm_1_V, Psi_h_1_v ] = Hao_3d_Partial_Derivatives_DSR_Extra(a,h,Tm,Sm,Sm_1_Mx,Sm_1_My,Sm_1_Tm,Sm_2_MxMx,Sm_2_MyMy,Sm_2_TmTm,Sm_2_MxMy,Sm_2_MxTm,Sm_2_MyTm,Sm_1_V,Sm_2_vTm,Psi_M,Psi_h)
%Derive the Partial derivatives for the Double Square Root Equation.
% ---input: ----
% a, h, psi_M: 2x1 vector
% Tm : scalar
% Sm, and related derivatives: 2x2 matrix
% ---Output: ----
% objective function: Fc  : 2x1 vector
% objective function's derivatives: Fc_1: 2x2 matrix
%Coded by Hao, Dec 1. 2016
%Modified by Hao to adapt to the 3D tomo scheme, Oct 5,2017
%
%%
% define the DSR equation
Ts = sqrt( 1/4*Tm^2 + (a-h).' * Sm * (a-h)); %sm: 2x2, a,h: 2x1; Ts: 1x1
Tr = sqrt( 1/4*Tm^2 + (a+h).' * Sm * (a+h)); %sm: 2x2, a,h: 2x1; Tr: 1x1
Td = Tr + Ts;

%% define the first order partial derivatives
%--- c1 ---
Ts_1_h = -1/Ts * Sm * (a-h); % Ts_1_h: 2x1
Tr_1_h =  1/Tr * Sm * (a+h); % Tr_1_h: 2x1

%--- c2 ---
Ts_1_a =  1/Ts * Sm * (a-h); % Ts_1_a: 2x1
Tr_1_a =  1/Tr * Sm * (a+h); % Tr_1_a: 2x1

%--- c3 ---
Ts_1_Mx =  1/(2*Ts) * (a-h).' * Sm_1_Mx * (a-h); % Ts_1_Mx: 1x1
Ts_1_My =  1/(2*Ts) * (a-h).' * Sm_1_My * (a-h); % Ts_1_My: 1x1
Ts_1_M  = [Ts_1_Mx; Ts_1_My];  % Ts_1_M: 2x1

Tr_1_Mx =  1/(2*Tr) * (a+h).' * Sm_1_Mx * (a+h); % Tr_1_Mx: 1x1
Tr_1_My =  1/(2*Tr) * (a+h).' * Sm_1_My * (a+h); % Tr_1_My: 1x1
Tr_1_M  = [Tr_1_Mx; Tr_1_My];  % Tr_1_M: 2x1

%--- c4 ---
Ts_1_Tm =  1/(Ts) * (Tm/4 + 1/2* (a-h).' * Sm_1_Tm * (a-h)); % Ts_1_Tm: 1x1 
Tr_1_Tm =  1/(Tr) * (Tm/4 + 1/2* (a+h).' * Sm_1_Tm * (a+h)); % Tr_1_Tm: 1x1

%% define the second order partial derivatives
%--- c5 ---
Ts_2_hh =  1/Ts * (Sm - Ts_1_h * Ts_1_h.' );  % Ts_2_hh: 2x2
Tr_2_hh =  1/Tr * (Sm - Tr_1_h * Tr_1_h.' );  % Tr_2_hh: 2x2

%--- c6 ---
Ts_2_aa =  1/Ts * (Sm - Ts_1_a * Ts_1_a.' );  % Ts_2_aa: 2x2
Tr_2_aa =  1/Tr * (Sm - Tr_1_a * Tr_1_a.' );  % Tr_2_aa: 2x2

%--- c7 ---
Ts_2_MxMx =  1/Ts * ( 1/2 * (a-h).' * Sm_2_MxMx * (a-h)  - Ts_1_Mx * Ts_1_Mx );  % Ts_2_MxMx: 1x1
Ts_2_MxMy =  1/Ts * ( 1/2 * (a-h).' * Sm_2_MxMy * (a-h)  - Ts_1_Mx * Ts_1_My );  % Ts_2_MxMy: 1x1
Ts_2_MyMy =  1/Ts * ( 1/2 * (a-h).' * Sm_2_MyMy * (a-h)  - Ts_1_My * Ts_1_My );  % Ts_2_MyMy: 1x1
Ts_2_MM   = [Ts_2_MxMx Ts_2_MxMy; Ts_2_MxMy Ts_2_MyMy]; % Ts_2_MM: 2x2

Tr_2_MxMx =  1/Tr * ( 1/2 * (a+h).' * Sm_2_MxMx * (a+h) - Tr_1_Mx * Tr_1_Mx );  % Tr_2_MxMx: 1x1
Tr_2_MxMy =  1/Tr * ( 1/2 * (a+h).' * Sm_2_MxMy * (a+h) - Tr_1_Mx * Tr_1_My );  % Tr_2_MxMy: 1x1
Tr_2_MyMy =  1/Tr * ( 1/2 * (a+h).' * Sm_2_MyMy * (a+h) - Tr_1_My * Tr_1_My );  % Tr_2_MyMy: 1x1
Tr_2_MM   = [Tr_2_MxMx Tr_2_MxMy; Tr_2_MxMy Tr_2_MyMy]; % Tr_2_MM: 2x2

%--- c8 ---
Ts_2_ha =  -1/Ts * (Sm + Ts_1_h * Ts_1_a.' ); % Ts_2_ha: 2x2
Tr_2_ha =   1/Tr * (Sm - Tr_1_h * Tr_1_a.' ); % Tr_2_ha: 2x2

%--- c9 ---
Ts_2_hMx =  -1/Ts * (Sm_1_Mx * (a-h) + Ts_1_Mx * Ts_1_h);  % Ts_2_hMx: 2x1
Ts_2_hMy =  -1/Ts * (Sm_1_My * (a-h) + Ts_1_My * Ts_1_h);  % Tr_2_hMy: 2x1
% Ts_2_hM  =  cat(2,Ts_2_hMx, Ts_2_hMy); % Ts_2_hM: 2x2 
Ts_2_hM  =  [ Ts_2_hMx.' ; Ts_2_hMy.']; % Ts_2_hM: 2x2

Tr_2_hMx =   1/Tr * (Sm_1_Mx * (a+h) - Tr_1_Mx * Tr_1_h ); % Sm_1_Mx: 2x2, Tr: 1x1, Ts_1_h: 2x1;  Tr_2_hMx: 2x1
Tr_2_hMy =   1/Tr * (Sm_1_My * (a+h) - Tr_1_My * Tr_1_h );
% Tr_2_hM  =   cat(2,Tr_2_hMx, Tr_2_hMy); % Tr_2_hM: 2x2
Tr_2_hM  =  [ Tr_2_hMx.' ; Tr_2_hMy.']; % Tr_2_hM: 2x2

%--- c10 ---
Ts_2_aMx =   1/Ts * (Sm_1_Mx * (a-h) - Ts_1_Mx * Ts_1_a ); % Ts_2_aMx: 2x1
Ts_2_aMy =   1/Ts * (Sm_1_My * (a-h) - Ts_1_My * Ts_1_a ); % Ts_2_aMy: 2x1
% Ts_2_aM  =   cat(2,Ts_2_aMx, Ts_2_aMy); % Ts_2_aM: 2x2
Ts_2_aM  =   [Ts_2_aMx.'; Ts_2_aMy.']; % Ts_2_aM: 2x2

Tr_2_aMx =   1/Tr * (Sm_1_Mx * (a+h) - Tr_1_Mx * Tr_1_a ); % Sm_1_Mx: 2x2, Ts_1_Mx:2x2,  Tr: 1x1, Tr_1_a: 2x1;  Tr_2_aMx: 2x1
Tr_2_aMy =   1/Tr * (Sm_1_My * (a+h) - Tr_1_My * Tr_1_a );
% Tr_2_aM  =   cat(2,Tr_2_aMx, Tr_2_aMy); % Tr_2_aM: 2x2
Tr_2_aM  =   [Tr_2_aMx.'; Tr_2_aMy.']; % Tr_2_aM: 2x2

%--- c11 ---
Ts_2_hTm =  -1/Ts * (Sm_1_Tm * (a-h) + Ts_1_Tm * Ts_1_h ); % Ts_2_hTm: 2x1
Tr_2_hTm =   1/Tr * (Sm_1_Tm * (a+h) - Tr_1_Tm * Tr_1_h ); % Tr_2_hTm: 2x1

%--- c12 ---
Ts_2_aTm =   1/Ts * (Sm_1_Tm * (a-h) - Ts_1_Tm * Ts_1_a ); % Ts_2_aTm: 2x1
Tr_2_aTm =   1/Tr * (Sm_1_Tm * (a+h) - Tr_1_Tm * Tr_1_a ); % Tr_2_aTm: 2x1

%--- c13 ---
Ts_2_MxTm =   1/Ts * (1/2 * (a-h).' * Sm_2_MxTm * (a-h) - Ts_1_Mx * Ts_1_Tm);  % Ts_2_MxTm: 1*1
Ts_2_MyTm =   1/Ts * (1/2 * (a-h).' * Sm_2_MyTm * (a-h) - Ts_1_Mx * Ts_1_Tm);  % Ts_2_MyTm: 1*1
Ts_2_MTm  =   [Ts_2_MxTm; Ts_2_MyTm]; % Ts_2_M: 2x1

Tr_2_MxTm =   1/Tr * (1/2 * (a+h).' * Sm_2_MxTm * (a+h) - Tr_1_Mx * Tr_1_Tm);  %  Tr_2_MxTm: 1*1
Tr_2_MyTm =   1/Tr * (1/2 * (a+h).' * Sm_2_MxTm * (a+h) - Tr_1_Mx * Tr_1_Tm);  %  Tr_2_MyTm: 1*1
Tr_2_MTm  =   [Tr_2_MxTm; Tr_2_MyTm]; % Tr_2_M: 2x1

%--- c14 ---
Ts_2_TmTm =  1/Ts * (1/4 + 1/2 * (a-h).' * Sm_2_TmTm * (a-h) - Ts_1_Tm * Ts_1_Tm); % Ts_2_TmTm: 1x1 
Tr_2_TmTm =  1/Tr * (1/4 + 1/2 * (a+h).' * Sm_2_TmTm * (a+h) - Tr_1_Tm * Tr_1_Tm); % Tr_2_TmTm: 1x1 


%% Derive a: Using Newton-Raphson Method to solving the (non)linear equations 

% define the function (consistance function)
Td_1_a  = Ts_1_a  + Tr_1_a;    % Td_1_a: 2x1 
Td_1_M  = Ts_1_M  + Tr_1_M;    % Td_1_M: 2x1
Td_1_Tm = Ts_1_Tm + Tr_1_Tm;   % Td_1_Tm 1x1

Fc = Td_1_a - Td_1_M - Td_1_Tm * Psi_M;  % Fc: 2x1

% define the first order derivative of this funciton
Td_2_aa  = Ts_2_aa  + Tr_2_aa;  % Td_2_aa:  2x2
Td_2_aM  = Ts_2_aM  + Tr_2_aM;  % Td_2_aM:  2x2
Td_2_aTm = Ts_2_aTm + Tr_2_aTm; % Td_2_aTm: 2x1 

Fc_1 = Td_2_aa - Td_2_aM - Psi_M * Td_2_aTm.'; % Fc: 2x2

%% Derive b: derive the partial derivative of dm/dv and dt/dv (Eq. 29-33)

%define the partial derivatives
Ts_1_V  = (a-h).' * Sm_1_V * (a-h) ./ (2*Ts) ;  % Ts_1_V:  1x1
Tr_1_V  = (a+h).' * Sm_1_V * (a+h) ./ (2*Tr) ;  % Tr_1_V:  1x1
Ts_2_va = Sm_1_V * (a-h)/Ts - (Ts_1_V) * Sm * (a-h)/(Ts^2); % Ts_2_va: 2x1
Tr_2_va = Sm_1_V * (a+h)/Tr - (Tr_1_V) * Sm * (a+h)/(Tr^2); % Tr_2_va: 2x1

Td_1_V  = Ts_1_V  + Tr_1_V;  % Td_1_V:   1x1
Td_2_va = Ts_2_va + Tr_2_va; % Td_2_va:  2x1

%define the euqations and solve dm/dv and dt/dv
A  = zeros(3,3);
b  = zeros(3,1);
Et = zeros(3,1);

A(1:2,1:2) = Td_2_aM - Td_2_aa;   % 2x2
A(1:2,3)   = Td_2_aTm;            % 2x1
A(3,1:2)   = (Td_1_M - Td_1_a).'; % 1x2  !! this has to be transposed !!
A(3,3)     = Td_1_Tm;             % 1x1

b(1:2,1)   = -Td_2_va;          % 2x1
b(3,1)     = -Td_1_V;           % 1x1

Et = A\b;
M_1_V  = Et(1:2,1); % 2x1
Tm_1_V = Et(3,1);   % 1x1

%% Derive c: derive the partial derivative of d_psih/dv (Eq. 35)

Td_2_hM   = Ts_2_hM   + Tr_2_hM;   % 2x2
Td_2_ha   = Ts_2_ha   + Tr_2_ha;   % 2x2
Td_2_MTm  = Ts_2_MTm  + Tr_2_MTm;  % 2x1
Td_2_hTm  = Ts_2_hTm  + Tr_2_hTm;  % 2x1
Td_2_TmTm = Ts_2_TmTm + Tr_2_TmTm; % 1x1

Ts_2_vh   = -(Sm_1_V * (a-h)/Ts  - Sm *(a-h)/(Ts)^2 * Ts_1_V);    % 2x1
Tr_2_vh   =  (Sm_1_V * (a+h)/Tr  - Sm *(a+h)/(Tr)^2 * Tr_1_V);    % 2x1
Ts_2_vTm  =  (a-h).' * Sm_1_V * (a-h)./2  *(-1/(Ts)^2)*Ts_1_Tm + (a-h).' * Sm_2_vTm * (a-h)./(2*Ts) ; % 1x1
Tr_2_vTm  =  (a+h).' * Sm_1_V * (a+h)./2  *(-1/(Tr)^2)*Tr_1_Tm + (a+h).' * Sm_2_vTm * (a+h)./(2*Tr) ; % 1x1

Td_2_vh   = Ts_2_vh  + Tr_2_vh;  % 2x1
Td_2_vTm  = Ts_2_vTm + Tr_2_vTm; % 1x1

%% derive the d_psih/dv (Eq. 35):

%select expression to derive d_Psi_h/dv, by default using plan_1
plan = 1;

if (plan==1)
    %(1)full expression
    if(size(Psi_h,1)>1)
        %in the case of half-offset vector(h) has two degree of freedom
        Psi_h_1_v = -(1/Td_1_Tm) * (Td_2_vh + Td_2_vTm * Psi_h + (Td_2_hM - Td_2_ha + (Td_2_MTm - Td_2_aTm)*(Psi_h.')) * M_1_V  + Tm_1_V * (Td_2_hTm + Td_2_TmTm * Psi_h)); % 2x1
    else
        %in the case of half-offset vector(h) has only one degree of freedom
        
%         %=== Old fomulation of d_h/d_sigma from Einar's old derivation (draft version: before Nov.01,2018)===
%         h_1_sigma = [h(1)/sqrt(h(1)^2 + h(2)^2);h(2)/sqrt(h(1)^2 + h(2)^2)]; % an unit vector: 2x1
        
        %=== New fomulation of d_h/d_sigma from Einar's new derivation (draft version: Nov.01,2018)===
%         n         = [1;0]; %assume the acqusion is along x direction, this has +to be changed in different case
        n         = [0;1]; % use this for Lundin's dataset
        h_1_sigma = (1./(n.' * h)).* h;
        
        Psi_sigma = Psi_h;
        Psi_h_1_v = -(1/Td_1_Tm) * (h_1_sigma.' * Td_2_vh  + Td_2_vTm * Psi_sigma + (h_1_sigma.' * (Td_2_hM - Td_2_ha) + (Td_2_MTm - Td_2_aTm).' * Psi_sigma) * M_1_V  + Tm_1_V * (h_1_sigma.' * Td_2_hTm + Td_2_TmTm * Psi_sigma)); % 1x1
    end
    
elseif (plan==2)
    %(2)modified expression, neglect both the dm_dv and dt_dv term to get faster convergence
    Psi_h_1_v = -(1/Td_1_Tm) * (Td_2_vh + Td_2_vTm * Psi_h);
    
else
    %(3)modified expression, neglect only the dm_dv term get faster convergence
    Psi_h_1_v = -(1/Td_1_Tm) * (Td_2_vh + Td_2_vTm * Psi_h + Tm_1_V * (Td_2_hTm + Td_2_TmTm * Psi_h));
    
end

end


