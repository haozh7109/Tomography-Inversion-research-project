function [ aper ] = Hao_Newton_Solver(H,Tm,Sm,Sm_1_M,Sm_1_Tm,Sm_2_MM,Sm_2_TmTm,Sm_2_MTm,Psi_m,diff_method)
%Using Newton-Raphson method to solving the (non)linear equation.
% based on Newton-Raphson algorithm to solve the unknow.

%Coded by Hao, Dec 2. 2016

%% define a intinital value of unknow
a0=0;

%% define the diffraction equation ( 1: dobule squre root, 2: single squre root,)  
method = diff_method;

%% define a stop condition (iterations or misfit error)
max_n  = 20;
misfit = 1e-5;

%% start the newton-raphson iteration

for in = 1:max_n
    
    %---start calculation---
    if (method ==1)
        [ f_a0, f_1_a0 ] =  Hao_Partial_Derivatives_DSR(a0,H,Tm,Sm,Sm_1_M,Sm_1_Tm,Sm_2_MM,Sm_2_TmTm,Sm_2_MTm,Psi_m);
    else
        [ f_a0, f_1_a0 ] =  Hao_Partial_Derivatives_SSR(a0,H,Tm,Sm,Sm_1_M,Sm_1_Tm,Sm_2_MM,Sm_2_TmTm,Sm_2_MTm,Psi_m);
    end
    
    %---update the unknow variable---
    a1 = a0 - f_a0/f_1_a0;
    
    if (method ==1)
        [ f_a1 ] = Hao_Partial_Derivatives_DSR(a1,H,Tm,Sm,Sm_1_M,Sm_1_Tm,Sm_2_MM,Sm_2_TmTm,Sm_2_MTm,Psi_m);
    else
        [ f_a1 ] = Hao_Partial_Derivatives_SSR(a1,H,Tm,Sm,Sm_1_M,Sm_1_Tm,Sm_2_MM,Sm_2_TmTm,Sm_2_MTm,Psi_m);
    end
    
    %---check if it satisfies misft---
    if(abs(f_a1) < misfit)
        a0 = a1;
        break
    elseif((abs(f_a1) > abs(f_a0)) && (abs(f_a0) < 10 * misfit))
        break
    end
    
    %---set x2 for the input of next iteration---
    a0 = a1;
    
end

%% output the searched aperture
aper = a0;

end

