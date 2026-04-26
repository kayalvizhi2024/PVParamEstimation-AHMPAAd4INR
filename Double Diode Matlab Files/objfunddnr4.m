function Fit=objfunddnr4(Xnew,C_Iter)
% objfunddnr4 - Fitness evaluation using 4th-order Newton-Raphson method
%
% Inputs:
%   Xnew    - Candidate solution vector [Rs, Rp, Iph, Io1, a1]
%   C_Iter  - Current iteration number
%
% Output:
%   Fit     - Root Mean Square Error (RMSE) between estimated and measured current
%
% This function estimates the output current of a PV cell using the
% Single Diode Model and refines it using a fourth-order Newton-Raphson method.
%==========================================================================

global dataset MaxIteration
Tc = 3.0615e+02;
q = 1.60217646e-19;
k = 1.3806503e-23;

a1=Xnew(6);              %Specified the diode ideality factor value from population
a2=Xnew(7);              %Specified the diode ideality factor value from population
Rs=Xnew(4);             %Specified the PV series resistance value from population
Rp=Xnew(5);             %Specified the PV parallel resistance value from population
Iph=Xnew(1);            %Specified the photo current value from population
Io1=Xnew(2);             %Specified the diode saturation current value from population
Io2=Xnew(3);             %Specified the diode saturation current value from population
             



dataset=[...
-0.2057	0.764
-0.1291	0.762
-0.0588	0.7605
0.0057	0.7605
0.0646	0.76
0.1185	0.759
0.1678	0.757
0.2132	0.757
0.2545	0.7555
0.2924	0.754
0.3269	0.7505
0.3585	0.7465
0.3873	0.7385
0.4137	0.728
0.4373	0.7065
0.459	0.6755
0.4784	0.632
0.496	0.573
0.5119	0.499
0.5265	0.413
0.5398	0.3165
0.5521	0.212
0.5633	0.1035
0.5736	-0.01
0.5833	-0.123
0.59	-0.21
];
VT=(k*Tc)/q;
Vp = dataset(:,1);            % Voltage points
Ie = dataset(:,2);            % Measured current points
N = length(Vp);               % Number of data points

% Adaptive damping factor
M = exp(-(5 * C_Iter / MaxIteration)^2.5);

F=Iph-Io1.*[exp((Vp+Ie.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ie.*Rs)./(a2*VT))-1]-((Vp+Ie.*Rs)./Rp)-Ie;
fd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ie.*Rs)./(a1*VT))-1))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ie.*Rs)./(a2*VT))-1))-(Rs/Rp)-1;
Ip=Ie-M.*(F./fd); 

FF=Iph-Io1.*[exp((Vp+Ip.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ip.*Rs)./(a2*VT))-1]-((Vp+Ip.*Rs)./Rp)-Ip;
fdd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ip.*Rs)./(a1*VT))))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ip.*Rs)./(a2*VT))-1))-(Rs/Rp)-1;
Ip=Ip-M.*(FF./fdd); 

FFF=Iph-Io1.*[exp((Vp+Ip.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ip.*Rs)./(a2*VT))-1]-((Vp+Ip.*Rs)./Rp)-Ip;
fddd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ip.*Rs)./(a1*VT))-1))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ip.*Rs)./(a2*VT))))-(Rs/Rp)-1;
Ip=Ip-M.*(FFF./fddd);

FFFF=Iph-Io1.*[exp((Vp+Ip.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ip.*Rs)./(a2*VT))-1]-((Vp+Ip.*Rs)./Rp)-Ip;
fdddd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ip.*Rs)./(a1*VT))))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ip.*Rs)./(a2*VT))-1))-(Rs/Rp)-1;

for i=1:length(Vp)
Ip(i)=Ie(i)-M*(F(i)/fd(i))-(F(i)^4/fdddd(i));
end

% Compute RMSE as fitness value
Fit = sqrt((1 / N) * sum((Ie - Ip).^2));

end
