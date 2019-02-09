clear all
clc

Name: Ysabelle Lam;
PID: A13868700;

%Assumptions%

% 1. The wing has a constant chord (rectangular wing)
% 2. Drag caused by wing only
% 3. Lift caused by wing and tail only (Ltot = Lwing + Ltail)
% 4. Aircraft is in steady level flight (Ltot = nW, T=D)
% 5. Tail is a symmetric airfoil (Cm,tail = 0)

%Data%
Data = xlsread('SE160A_Project1_InputFile.xlsx','Input');

%Aircraft Geometry Definition%

wingstation = Data(1,1);
tailstation = Data(2,1);
ww = Data(3,1);
S = Data(4,1);
wingspan = Data(5,1);
c = Data(6,1);
St = Data(7,1);
bt = Data(8,1);
ct = Data(9,1);

%Aircraft Aerodynamic Definintion%

lift_slope = Data(14,1);
ao = Data(15,1);
CLmaxpos = Data(16,1);
stallanglepos = Data(17,1);
CLmaxneg = Data(18,1);
stallangleneg = Data(19,1);
Cmo = Data(20,1);
Cdo = Data(21,1);
K = Data(22,1);
do = Data(23,1);
d_10 = Data(24,1);

%Aircraft Performance%

V_stall = Data(29,1);
V_stall = V_stall*1.46667; %mph to ft/sec
V_cruise = Data(30,1);
V_cruise = V_cruise*1.46667; %mph to ft/sec
V_dive = Data(31,1);
V_dive = V_dive*1.46667; %mph to ft/sec
npos = Data(32,1);
nneg = Data(33,1);

%Aircraft Weigth and Balance%

w_aircraft = Data(38,1);
l_aircraft = Data(38,2);
w_fuel = Data(39,1);
l_fuel = Data(39,2);
w_pilot = Data(40,1);
l_pilot = Data(40,2);
w_copilot = Data(41,1);
l_copilot = Data(41,2);
w_pass = Data(42,1);
l_pass = Data(42,2);
w_luggage = Data(43,1);
l_luggage = Data(43,2);

%Aircraft Analysis Study%

density = Data(48,1); 
analysis = Data(49,1);
max_weight = Data(50,1);
min_CG = Data(51,1);
max_CG = Data(52,1);

%Calculations%

%Weight & Center of Gravity Calculations%

w_total = w_aircraft + w_copilot + w_fuel + w_luggage + w_pass + w_pilot; 
w_l = (w_aircraft*l_aircraft) + (w_copilot*l_copilot) + ...
    (w_fuel*l_fuel) + (w_luggage*l_luggage) + (w_pass*l_pass) + ...
    (w_pilot*l_pilot);
center_gravity = w_l/w_total;

if w_total > max_weight
    f = msgbox('WARNING: Aircraft weight exceeds allowable limits', ...
        'Warning');
end

if center_gravity > max_CG
    f = msgbox('WARNING: Center of Gravity is outside allowable limits', ...
        'Warning');
else if center_gravity < min_CG
    f = msgbox('WARNING: Center of Gravity is outside allowable limits', ...
        'Warning');
    end
end

%Aerodynamic Calculations%

V_PHAA = V_stall * sqrt(npos);

V_NHAA = V_stall * sqrt(-nneg);

V_PLAA = V_dive;

V_NLLA = V_dive;

n_pos_dive = 1 + (.5*(density)*(S/w_total)*(lift_slope*57.2957549575)*(25)*(V_dive));

n_pos_cruise = 1 + (.5*(density)*(S/w_total)*(lift_slope*57.2957549575)*(50)*(V_cruise));

n_neg_dive = 1 + (.5*(density)*(S/w_total)*(lift_slope*57.2957549575)*(-25)*(V_dive));

n_neg_cruise = 1 + (.5*(density)*(S/w_total)*(lift_slope*57.2957549575)*(-50)*(V_cruise));


%Airspeed is dependent on input data

if analysis == 1
    Airspeed = V_PHAA;
else if analysis == 2
        Airspeed = V_PLAA;
else if analysis == 3
        Airspeed = V_NHAA;
else if analysis == 4
        Airspeed = V_NLAA;
else if analysis == 5
        Airspeed = V_cruise;
else if analysis == 6
        Airspeed = V_dive;
else if analysis == 7
        Airspeed = V_cruise;
    else analysis == 8
        Airspeed = V_dive;
    end
    end
    end
    end
    end
    end
    end


%LoadFactor is dependent on input data
if analysis == 1
    n = npos;
else if analysis == 2
    n = npos;
else if analysis == 3
    n = nneg;
else if analysis == 4
    n = neg;
else if analysis == 5
    n = n_pos_cruise;
else if analysis == 6
    n = n_pos_dive;
else if analysis == 7
    n = n_neg_cruise;
    else analysis == 8
    n = n_neg_dive;
    end
    end
    end
    end
    end
    end
end

%AerodynamicMoment%

%Fx = T - D
%Fy = 0 = L + Lt - Wn
%M = 0 = Ma - Lt(tail station - wing station in feet)

q = .5*density*(Airspeed^2);
Moment = q*S*c*Cmo;

%Lift
fy_matrix = [1,1;((wingstation-center_gravity)/12),((tailstation-center_gravity)/12)];
moment_matrix = [n*w_total;Moment];
lift_matrix = inv(fy_matrix)*moment_matrix;

%Tail Lift
TailLift = lift_matrix(1);

%Wing Lift
WingLift = lift_matrix(2);

CL = (2*n*w_total)/(density*(Airspeed^2)*(S));

%Drag%
Cd = Cdo + K*(CL^2);
Drag = Cd*q*S;

%Thrust%
% Here we assumed that drag = thrust%
Thrust = Drag;

%AngleofAttack%
angle_attack = CL/lift_slope + ao;

%V-N Diagram%

% Axis

x2 =  linspace(0,300);
y = 0 + (0*x2);

y0 = linspace(-4,7);
x0 = 0 + (0*y0);
hold;

plot(x2,y,'k')

%V-n diagram main graph

a = linspace(V_PHAA,V_dive);
y1 = npos + (0*a);
p1 = plot(a,y1,'k')

b = linspace(V_stall,V_cruise);
y2 = nneg + (0*b);
plot(b,y2,'k')

b1 = linspace(0,npos);
x1 = V_dive + (0*b1);
plot(x1,b1,'k') 

g = linspace(0,V_stall);
y9 = (npos/V_PHAA^2)*(g).^2;
plot(g,y9,'k')

h = linspace(0,V_stall);
y10 =(nneg/(V_stall^2))*(h).^2;  %nneg/V_stall%
plot(h,y10,'k')

g1 = linspace(V_stall,V_PHAA);
y91 = (npos/V_PHAA^2)*(g1).^2;
plot(g1,y91,'k')

t = linspace(V_cruise,V_dive);
y11 = ((0-nneg)/(V_dive-V_cruise))*t - (((0-nneg)*V_dive)/(V_dive-V_cruise));
plot(t,y11,'k')

%Gust Curves

d = linspace(0,V_cruise);
y3 = 1 + (((n_pos_cruise)/V_cruise)*d);
p2 = plot(d,y3,'--r');

e = linspace(V_cruise,V_dive );
y4 = ((n_pos_cruise) - ((n_pos_dive-n_pos_cruise)/(V_dive-V_cruise))*V_cruise...
    +(((n_pos_dive-n_pos_cruise)/(V_dive-V_cruise))*e))+1;
plot(e,y4,'--r')

e1 = linspace(0,V_dive);
y5 = 1 + (((n_pos_dive)/V_dive)*e1);
plot(e1,y5,'--r')

y6 = 1 + (((n_neg_cruise)/V_cruise) * d);
p3 = plot(d,y6,'--r'); 

y7 = ((((n_neg_dive-n_neg_cruise)/(V_dive-V_cruise)) * e)- ...
    (((n_neg_dive-n_neg_cruise)/(V_dive-V_cruise))*V_cruise)...
    + n_neg_cruise)+1;
plot(e,y7,'--r')

y8 = 1 + (((n_neg_dive)/V_dive) * e1);
plot(e1,y8,'--r')

z = linspace(y4(end),y7(end));
x1 = V_dive + (0*z);
plot(x1,z,'--r') 

title('V-n Diagram')
xlabel('Velocity (mph)')
ylabel('Load Factor (n)')
legend({'V-n Diagram','Gust Envelope'},'Location','northwest')

%add color%

%Distributions 

%Taylor series of lift distribution
%alpha = 1 in the case of a rectangular wing

btay1 = (-wingspan/2);
btay2 = (wingspan/2);

x = linspace(btay1,btay2);

distributed_lift = (WingLift/wingspan)*(((2/pi)+.5)-((2*x./wingspan).^2).*(1/pi))-(((2*x./wingspan)).^4).*(1/4*pi)-...
    ((2*x./wingspan).^6).*(1/8*pi)-((2*x./wingspan).^8).*(5/64*pi)-((2*x./wingspan).^10).*(7/128*pi);


%Wing Moment
dw = ones(100);
distributed_moment = (c*.5*density*(Airspeed^2)*Cmo).*dw;


%distributed_drag
fun = d_10.*(2*(x.^2)./wingspan).^10;
fun = fun + do;
coeff = (.5*density*(Airspeed^2)*do*(1+(-2/wingspan))*(Cdo+K*(CL^2)));
distributed_drag = fun.*coeff;

%distributed weight
distributed_weight = (n*w_total/wingspan).*dw;

%graphs
figure
subplot(2,2,1);
plot(x, distributed_lift)
title('Lift Distribution')
xlabel('Distance Along Span(ft)')
ylabel('Lift')

subplot(2,2,2);
plot(x, distributed_moment)
title('Moment Distribution')
xlabel('Distance Along Span(ft)')
ylabel('Moment')

subplot(2,2,3);
plot(x, distributed_drag)
title('Drag Distribution')
xlabel('Distance Along Span(ft)')
ylabel('Drag')

subplot(2,2,4);
plot(x, distributed_weight)
title('Weight Distribution')
xlabel('Distance Along Span(ft)')
ylabel('Weight')

%For grading
Results = [w_total, center_gravity, Airspeed, n, Moment, WingLift, TailLift, ...
    Drag, Thrust, angle_attack ];