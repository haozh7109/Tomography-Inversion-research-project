%% Job used to fitting the derived finite Surfaces to the seismic
%==========================================================================
% Job used to tuning the derived finite Surfaces with the seismic gathers
% coded by Hao, 29th Jan. 2018

%% load the input surface(curve) and CMP gather
hor_input   = load('tomo3d_test_calculated_horizon.mat');
hor_t_cell  = hor_input.hor_t;
hor_t0_cell = hor_input.hor_t0;
hor_x_cell  = hor_input.hor_offset;
cmp_input   = load('tomo3d_test_cmp.mat');
cmp         = cmp_input.cmp_gather;

QC_hor      = 1:5; %define the QC horizon(s)

%% plot the input surface
figure;
% plot the horizons
for ihor = 1:5
    subplot(1,2,1)
    plot(hor_t_cell{ihor})
    hold on
end
set(gca,'Ydir','reverse')
title('simulated finite offset horizons')

% plot the cmp gather and horizons
subplot(1,2,2)
plotseis(cmp(1:2000,:),1:2000,1:94,1,[5 10])
for ihor = 1:5
    hold on
    subplot(1,2,2)
    scatter(1:94,round(hor_t_cell{ihor}./0.002),30,'o');
end
title('horizons overlaied with cmp gather')

%% (1) fit the surface with 4th order reflection travel time equation (Alkhalifah and Tsvankin 1995)

% define the cell to save the fitting parameters
fit1_func     = cell(1,5);
fit1_para_v   = cell(1,5);
fit1_para_eta = cell(1,5);


figure;
for ihor = QC_hor%1:5
    
    % assign the input variables
    x  = hor_x_cell{ihor};
    t0 = hor_t0_cell{ihor};
    t  = hor_input.hor_t{ihor};
    
    % -- Fitting the NMO velocity and Eta parameters on the one go ---
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[1300,0],...
        'Upper',[6000,Inf],...
        'StartPoint',[1 1]);
    ft = fittype('sqrt(t0^2 + x.^2/v^2 - 2*eta*x.^4 ./( v^2 * (t0^2 * v^2 + (1+2*eta).* x.^2)))','independent',{'x'},'dependent',{'t'},'coefficients',{'v','eta'},'problem',{'t0'},'options',fo);% 4th order
    
    % fit the curve and make the plot
    [curve1,gof1] = fit(x',t',ft,'problem',{t0});
    plot(curve1,x,t)
    hold on;
    
%     % -- 1st fitting to derive the NMO velocity ---
%     fo = fitoptions('Method','NonlinearLeastSquares',...
%         'Lower',[0],...
%         'Upper',[Inf,Inf],...
%         'StartPoint',[1000]);
%     ft = fittype('sqrt(t0^2 + x^2/v^2)','independent',{'x'},'dependent',{'t'},'coefficients',{'v'},'problem',{'t0'},'options',fo);% 2nd order
%     [curve2,gof2] = fit(x',t',ft,'problem',{t0});
%     figure;plot(curve2,x,t)
%     
%     % -- 2nd fitting to derive the Eta parameter ---
%     
%     fo = fitoptions('Method','NonlinearLeastSquares',...
%         'Lower',[-Inf],...
%         'Upper',[Inf],...
%         'StartPoint',[0]);
%     
%     ft = fittype('sqrt(t0^2 + x^2/v^2 - 2*eta*x^4/(v^2*(t0^2*v^2+(1+2*eta)*x^2)))','independent',{'x'},'dependent',{'t'},'coefficients',{'eta'},'problem',{'t0','v'},'options',fo);% 4th order
%     [curve3,gof3] = fit(x',t',ft,'problem',{t0,curve2.v});
%     figure;plot(curve3,x,t)

    % save the fit parameters to output 
    fit1_func{ihor}     = curve1;
    fit1_para_v{ihor}   = curve1.v;
    fit1_para_eta{ihor} = curve1.eta;
    
    if(ihor ==1)
    fprintf('============ 1st pass fitting ===============\n')
    end
    fprintf('=========== fitting derived parameters for horizon %d: V = %f m/s, eta = %f ==============\n',ihor, curve1.v, curve1.eta);
end
set(gca,'Ydir','reverse')
title(' ====== finite offset horizons and fitted curve ===== ')


%% (2) perturbation of the fitted surface to further match with seismic events using semblance approach

% define the cell to save the fitting parameters
fit2_para_v    = cell(1,5);
fit2_para_eta  = cell(1,5);
fit2_semblance = cell(1,5);
hor_tnew_cell  = cell(1,5);

for ihor = QC_hor%1:5
    
    % assign the input variables
    x  = hor_x_cell{ihor};
    t0 = hor_t0_cell{ihor};
    
    % define the searching parameters of velocity and eta
    v_min    = fit1_para_v{ihor}*(1-0.5);
    v_max    = fit1_para_v{ihor}*(1+0.5);
    n_v      = 100;
    
    eta_min  = fit1_para_eta{ihor}*(1-0.5);
    eta_max  = fit1_para_eta{ihor}*(1+0.5);
    n_eta    = 100;
    
    v_test   = linspace(v_min,v_max,n_v);
    eta_test = linspace(eta_min,eta_max,n_eta);
   
    % define the searching need parameters
    dt       = 0.002;                            % sample ratio of seismic
    nt       = size(cmp,1);                      % sample length of seismic
    L        = 10;                               % half vertical length of summing window
    nh       = length(x);                        % get the number of offsets
    H_win    = hanning(2*L+1)*ones(1,length(x)); %define a hanning window(size:(2L+1)*Noff)for summing
    Semb_mat = zeros(n_v,n_eta);
    
    cmp_flag_mat = ones(size(cmp));              %define a flag mat to differentiate 'zero' and 'non-zero' values
    cmp_flag_mat(abs(cmp)<1e-10) = 0;
    
    % start the searching
    for iv = 1:n_v
        for ieta = 1:n_eta
            %get current veloity and eta
            v   = v_test(iv);
            eta = eta_test(ieta);
            
            % derive the central travel time
            time = sqrt(t0^2 + x.^2/v^2 - 2*eta*x.^4 ./( v^2 * (t0^2 * v^2 + (1+2*eta).* x.^2)));
            
            % derive the summing window for every trace
            win_stat = round(time./dt) - L + 1;
            win_end  = round(time./dt) + L + 1;
            amp_mat  = zeros(2*L+1,nh);
            flag_mat = zeros(2*L+1,nh);
            
            %extract the amplitude from the data
            for ih = 1:nh
                if (win_stat(ih) >= 1 && win_end(ih) <= nt)
                    amp_mat(:,ih) = cmp(win_stat(ih):win_end(ih),ih);
                    flag_mat(:,ih)= cmp_flag_mat(win_stat(ih):win_end(ih),ih);  
                end
            end
            
            % apply the hanning window and calculate the semblance value
            amp_mat           = amp_mat.* H_win;
            s1                = sum((sum(amp_mat,2)).^2);
            s2                = sum(sum((amp_mat).^2),2);
%             Semb_mat(iv,ieta) = (1/nh)*(s1/s2);
            
            % weight the semblance by totoal number of non-zero values
            n_nan_zero = sum(cmp_flag_mat(:));
            Semb_mat(iv,ieta) = (1/n_nan_zero)*(s1/s2);
        end
    end
    
    % save the fit parameters to output
    fit2_semblance{ihor} = Semb_mat;

    % find the max semblance
    [max_val,abs_ind] = max(Semb_mat(:));
    [vel_row,eta_col] = ind2sub(size(Semb_mat),abs_ind);
    
    % save the new derived velocity and eta
    fit2_para_v{ihor}    = v_min    + (vel_row-1) .* (  v_max-v_min)/n_v;
    fit2_para_eta{ihor}  = eta_min  + (eta_col-1) .* (eta_max-eta_min)/n_eta;
    
    % save the new curves
    v_new   = fit2_para_v{ihor};
    eta_new = fit2_para_eta{ihor};
    t_new   = sqrt(t0^2 + x.^2/v_new^2 - 2*eta_new*x.^4 ./( v_new^2 * (t0^2 * v_new^2 + (1+2*eta_new).* x.^2)));
    hor_tnew_cell{ihor}  = t_new;
    
    % update the slope_h based on the new fitted curves 
    % -- option-1, use symblic expression to derive partial derivatives
    syms X
    T       = symfun(sqrt(t0^2 + X.^2/v_new^2 - 2*eta_new * X.^4 ./( v_new^2 * (t0^2 * v_new^2 + (1+2*eta_new).* X.^2))),[X]);
    T_X     = diff(T);
    t_x_new = double(T_X(x)); 
    
    % -- option-2, use curve fitting to derive partial derivatives
    [curve2,gof2] = fit(x',t_new',ft,'problem',{t0});
    t_x_old       = differentiate(fit1_func{ihor},x');
    t_x_new2      = differentiate(curve2,x');
    
    %-- compare the old slope and current slope
    figure;plot(x,t_x_new,'r');hold on;plot(x,t_x_new2,'g');hold on;plot(x,t_x_old,'b');
    title(['slope comparison of surface:' num2str(ihor)]);set(gca,'Ydir','Reverse');
    
    % print the new fitting result
    if(ihor ==1)
    fprintf('============ 2nd pass fitting ===============\n')
    end
    fprintf('=========== fitting derived parameters for horizon %d: V = %f m/s, eta = %f ==============\n',ihor, v_new, eta_new);
  
end

%% QC the result
disp('===== Fitting completed ===== ')

figure;
% plot the cmp gather and original horizons
subplot(1,2,1)
plotseis(cmp(1:2000,:),1:2000,1:94,1,[5 10])
for ihor = QC_hor%1:5
    hold on
    subplot(1,2,1)
    scatter(1:94,round(hor_t_cell{ihor}./0.002),100,'*');
end
title('original horizons overlaied with cmp gather')

% plot the cmp gather and fitting updated horizons
subplot(1,2,2)
plotseis(cmp(1:2000,:),1:2000,1:94,1,[5 10])
for ihor = QC_hor%1:5
    hold on
    subplot(1,2,2)
    scatter(1:94,round(hor_tnew_cell{ihor}./0.002),100,'*');
end
title('New derived horizons overlaied with cmp gather')


















