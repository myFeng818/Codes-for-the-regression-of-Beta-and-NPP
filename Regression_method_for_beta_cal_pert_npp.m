clc;
clear all;
close all;


%% To calculate the feedback parameter beta = dNPP/dCO2
%  First, estimate the sensitivities of NPP to P and T, with detrended
%  values, i.e., NPPde ~ f(Tde, Pde)
% Next, remove the effect of T and P with the sensitivities estimated 
% with detrended values, i.e., Residual = NPP - PSens_de*P - TSens_de*T;
% Finally, estimate the sensitivities of Residual to CO2, i.e.,Residual ~ f(CO2)

%% 1.NPPde ~ f(Tde, Pde)

Model_list = {'CESM2',...
              'CESM2-WACCM',...
              'EC-Earth3-Veg',...
              'EC-Earth3-Veg-LR',...
              'MIROC-ES2L',...
              'MPI-ESM-1-2-HAM',...
              'MPI-ESM1-2-LR',...
              'NorESM2-LM',...
              'NorESM2-MM',...
              'UKESM1-0-LL'};

% Data load: Precipitation, Temperature, NPP and CO2 time series; from
% 1850-2014
Model_P_Boreal = xlsread('Precipitation_annual_mean_Boreal_USE.xlsx');
Model_T_Boreal = xlsread('Temperature_annual_mean_Boreal_USE.xlsx');
Model_NPP_Boreal = xlsread('NPP_annual_mean_Boreal.xlsx');
Model_CO2 = xlsread('CO2_annual_mean_1960-2014.xlsx');

% Change the Unit of P, T during 1960-2014
Model_P_Boreal_USE = Model_P_Boreal(111:end,2:end).*86400.*365;
Model_T_Boreal_USE = Model_T_Boreal(111:end,2:end)-273.5;
% Calculate the percentage increase of NPP relative to the mean NPP during 1950-1959
Model_NPP_Boreal_USE = (Model_NPP_Boreal(111:end,2:end)  - repmat(mean(Model_NPP_Boreal(91:110,2:end),1),55,1))./repmat(mean(Model_NPP_Boreal(91:110,2:end),1),55,1);

Model_P_diff = mean(Model_P_Boreal_USE(end-9:end,:),1) - mean(Model_P_Boreal_USE(1:10,:),1);
Model_T_diff = mean(Model_T_Boreal_USE(end-9:end,:),1) - mean(Model_T_Boreal_USE(1:10,:),1); 
Model_NPP_diff = mean(Model_NPP_Boreal_USE(end-9:end,:),1) - mean(Model_NPP_Boreal_USE(1:10,:),1); 

sens_long_P = Model_NPP_diff./Model_P_diff;
sens_long_T = Model_NPP_diff./Model_T_diff;

for imodel =1:10
    % Detrend the P, T and NPP (%)
    Model_P_Boreal_de(:,imodel) = detrend(Model_P_Boreal_USE(:,imodel));
    Model_T_Boreal_de(:,imodel) = detrend(Model_T_Boreal_USE(:,imodel));
    Model_NPP_Boreal_de(:,imodel) = detrend(Model_NPP_Boreal_USE(:,imodel));
    
    % Regression of NPP with P and T, with detrended values
    Xmatrix = [ones(55,1) Model_P_Boreal_de(:,imodel) Model_T_Boreal_de(:,imodel)];
    Y = Model_NPP_Boreal_de(:,imodel);
    [b,bint,r,rint,stats]=regress(Y,Xmatrix);
    
    Psens_de(imodel,1) = b(2);
    Tsens_de(imodel,1) = b(3);
    R2_de(imodel,1) = stats(1);
    p_de(imodel,1) = stats(3);
    
end

%% 2.Residual ~ f(CO2)
clear b bint r rint stats
figure,
for imodel =1:10
    % The original (not detrended) series of P, T, CO2, and NPP
    Model_P_Boreal_org(:,imodel) = Model_P_Boreal_USE(:,imodel);
    Model_T_Boreal_org(:,imodel) = Model_T_Boreal_USE(:,imodel);
    Model_NPP_Boreal_org(:,imodel) = Model_NPP_Boreal_USE(:,imodel);
    Model_CO2_org(:,imodel) = Model_CO2(:,2);
    
    % Remove the change caused by P and T using sensitivities from
    % detrended regression analysis
    Residual = Model_NPP_Boreal_org(:,imodel) - Psens_de(imodel,1).*Model_P_Boreal_org(:,imodel) - Tsens_de(imodel,1).*Model_T_Boreal_org(:,imodel);
    
    % Regression of Residual with CO2
    Xmatrix = [ones(55,1) Model_CO2_org(:,imodel)];
    Y = Residual;
    [b,bint,r,rint,stats]=regress(Y,Xmatrix);

    beta_new(imodel,1) = b(2);
    R2_new(imodel,1) = stats(1);
    p_new(imodel,1) = stats(3);
    
    % Simple plot to show the regressed relationship
    subplot(4,4,imodel),
    
    scatter(Model_CO2_org(:,imodel), Residual,'filled'); hold on;
    x1fit = min(Model_CO2_org(:,imodel)):0.1:max(Model_CO2_org(:,imodel));
    yfit = x1fit.*b(2) + b(1);
    plot(x1fit,yfit,'k-');

    xlabel('CO_2','Fontsize',5);
    ylabel('Residual','Fontsize',5);
    
    title(Model_list{imodel},'Fontsize',5);
     set(gca,'Fontsize',5);
%      view(50,10);
     hold off
    
end

% The calculated beta_NPP, response of plant growth to elevated CO2
% concentration
beta_matrix(:,1) = beta_new(:,1);

% The R2 and p values of the two regression models
statistics_index(:,1) = R2_de(:,1);
statistics_index(:,2) = R2_new(:,1);
statistics_index(:,3) = p_de(:,1);
statistics_index(:,4) = p_new(:,1);

% Save the beta, R2 and p values
xlswrite('Beta_Boreal_pert_npp.xlsx',beta_matrix);
xlswrite('R2_p_Boreal_pert_npp.xlsx',statistics_index);









