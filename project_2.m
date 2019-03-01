clear all
clc

%Data%
data = xlsread('SE160A_Project2_InputFile.xlsx','Stress');

sigma_xx = data(1);
sigma_yy = data(2);
sigma_zz = data(3);
tau_yz = data(4);
tau_xz = data(5);
tau_xy = data(6);

tension_y = data(11);
tension_u = data(12);
compression_y = data(13);
compression_u = data(14);
shear_y = data(15);
shear_u = data(16);
FSy = data(17);
FSu = data(18);

%%Calculations%%

%Allowable Stresses%
tension_allow = tension_y/FSy;
compression_allow = compression_y/FSy;
shear_allow = shear_y/FSy;


%Principle Stresses%
stress_matrix = [sigma_xx, tau_xy, tau_xz; tau_xy, sigma_yy, tau_yz;...
    tau_xz, tau_yz, sigma_zz];

e = eig(stress_matrix);

p_1 = max(e);
p_2 = median(e);
p_3 = min(e);


%Max Shear Stress%
shear_max = .5*(p_1-p_3);


%Mohr Circle%
angles=0:0.01:2*pi;

tau_1 = .5*(p_1-p_3);
tau_2 = .5*(p_1-p_2);
tau_3 = .5*(p_2-p_3);

c_1=(.5*(p_1+p_3));
c_2=(.5*(p_1+p_2));
c_3=(.5*(p_2+p_3));

circle_1=[c_1 + (tau_1*cos(angles')), tau_1*sin(angles')];
circle_2=[c_2 + (tau_2*cos(angles')), tau_2*sin(angles')];
circle_3=[c_3 + (tau_3*cos(angles')), tau_3*sin(angles')];

plot(circle_1(:,1),circle_1(:,2),circle_2(:,1),circle_2(:,2),...
    circle_3(:,1),circle_3(:,2));

axis 'equal';
patch(circle_1(:,1),circle_1(:,2),'cyan') ;
patch(circle_2(:,1),circle_2(:,2),'white') ;
patch(circle_3(:,1),circle_3(:,2),'white') ;

text(p_1,0,'\sigma_1','fontsize',10)
text(p_2,0,'\sigma_2','fontsize',10)
text(p_3,0,'\sigma_3','fontsize',10);

title('3D Mohrs Circle Stress States')
xlabel('\sigma(ksi)')
ylabel('\tau(ksi)')

legend({'\sigma_1 & \sigma_3 circle','\sigma_1 & \sigma_2 circle',...
    '\sigma_2 & \sigma_3 circle','Possible Stress States',...
    'Impossible Stress States', 'Impossible Stress States'},'Location',...
    'southoutside','Orientation','vertical')


%Margin of Safety%

MS_tresca = min([(tension_allow/p_1), (compression_allow/p_3), (shear_allow/shear_max)]); 
MS_tresca = abs(MS_tresca) - 1;

MS_rankine = min([(tension_allow/p_1), (compression_allow/p_3)]); 
MS_rankine = abs(MS_rankine) - 1;

sigma_eff = sqrt(.5*((p_1-p_3)^2 + (p_2-p_3)^2 + (p_1-p_2)^2));
MS_von_Mises = (tension_allow/sigma_eff) - 1;


%Max Stress States%

Max_Stress_State_Tresca = (MS_tresca +1).*stress_matrix;

Max_Stress_State_Rankine = (MS_rankine +1).*stress_matrix;

Max_Stress_State_von_Mises = (MS_von_Mises +1).*stress_matrix;


%Results%

Results = [tension_allow, compression_allow, shear_allow, p_1, p_2, p_3, shear_max,...
    MS_rankine, MS_tresca, MS_von_Mises];
