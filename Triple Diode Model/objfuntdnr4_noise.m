function Fit=objfuntdnr4_noise(Xnew,C_Iter)
global dataset MaxIteration

q = 1.60217646e-19;
k = 1.3806503e-23;
Tc = 3.0615e+02;

a1=Xnew(3);              %Specified the diode ideality factor value from population
a2=Xnew(4);              %Specified the diode ideality factor value from population
a3=Xnew(5);              %Specified the diode ideality factor value from population
Rs=Xnew(1);             %Specified the PV series resistance value from population
Rp=Xnew(2);             %Specified the PV parallel resistance value from population
Iph=Xnew(6);            %Specified the photo current value from population
Io1=Xnew(7);             %Specified the diode saturation current value from population
Io2=Xnew(8);             %Specified the diode saturation current value from population
Io3=Xnew(9);             %Specified the diode saturation current value from population


% 
dataset=[...
-0.2057	0.764	7.64E-01	7.64E-01	7.64E-01	7.58E-01	7.36E-01	8.11E-01
-0.1291	0.762	7.62E-01	7.62E-01	7.63E-01	7.67E-01	7.57E-01	8.03E-01
-0.0588	0.7605	7.61E-01	7.61E-01	7.60E-01	7.56E-01	7.28E-01	7.49E-01
0.0057	0.7605	7.61E-01	7.61E-01	7.60E-01	7.58E-01	7.34E-01	8.24E-01
0.0646	0.76	7.60E-01	7.60E-01	7.61E-01	7.65E-01	7.66E-01	8.17E-01
0.1185	0.759	7.59E-01	7.59E-01	7.60E-01	7.52E-01	7.66E-01	6.98E-01
0.1678	0.757	7.57E-01	7.57E-01	7.58E-01	7.59E-01	7.57E-01	7.19E-01
0.2132	0.757	7.57E-01	7.57E-01	7.57E-01	7.50E-01	7.20E-01	8.25E-01
0.2545	0.7555	7.56E-01	7.55E-01	7.55E-01	7.52E-01	7.66E-01	8.25E-01
0.2924	0.754	7.54E-01	7.54E-01	7.54E-01	7.55E-01	7.71E-01	7.99E-01
0.3269	0.7505	7.50E-01	7.51E-01	7.50E-01	7.47E-01	7.55E-01	7.98E-01
0.3585	0.7465	7.47E-01	7.47E-01	7.46E-01	7.44E-01	7.37E-01	7.02E-01
0.3873	0.7385	7.39E-01	7.39E-01	7.39E-01	7.34E-01	7.57E-01	7.25E-01
0.4137	0.728	7.28E-01	7.28E-01	7.28E-01	7.31E-01	7.61E-01	7.59E-01
0.4373	0.7065	7.07E-01	7.07E-01	7.06E-01	7.07E-01	7.13E-01	6.39E-01
0.459	0.6755	6.76E-01	6.76E-01	6.76E-01	6.72E-01	6.94E-01	6.09E-01
0.4784	0.632	6.32E-01	6.32E-01	6.32E-01	6.29E-01	6.39E-01	5.93E-01
0.496	0.573	5.73E-01	5.73E-01	5.73E-01	5.71E-01	5.52E-01	6.08E-01
0.5119	0.499	4.99E-01	4.99E-01	4.99E-01	5.04E-01	4.93E-01	4.52E-01
0.5265	0.413	4.13E-01	4.13E-01	4.13E-01	4.15E-01	4.04E-01	3.93E-01
0.5398	0.3165	3.17E-01	3.16E-01	3.16E-01	3.16E-01	3.31E-01	2.95E-01
0.5521	0.212	2.12E-01	2.12E-01	2.12E-01	2.13E-01	2.08E-01	2.19E-01
0.5633	0.1035	1.04E-01	1.04E-01	1.04E-01	1.03E-01	1.03E-01	1.02E-01
0.5736	-0.01	-1.00E-02	-1.00E-02	-1.00E-02	-1.01E-02	-9.85E-03	-9.27E-03
0.5833	-0.123	-1.23E-01	-1.23E-01	-1.23E-01	-1.22E-01	-1.24E-01	-1.20E-01
0.59	-0.21	-2.10E-01	-2.10E-01	-2.10E-01	-2.12E-01	-2.11E-01	-2.19E-01
];

VT=(k*Tc)/q;
Vp=dataset(:,1);
Ie=dataset(:,6); %change the column number to simulate the effect of noise

Ip=zeros(size(Vp));
N = length(Vp); % define Rp, Rs, Vth, a – Ideality factor, Io, Ipv, … 

%% 4th order NR

M=exp(-(5*C_Iter/MaxIteration)^2.5);
% M1=1-(C_Iter/(MaxIteration));
% M2=1-(C_Iter/(MaxIteration));
% M=1;
F=Iph-Io1.*[exp((Vp+Ie.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ie.*Rs)./(a2*VT))-1]...
    -Io3.*[exp((Vp+Ie.*Rs)./(a3*VT))-1]-((Vp+Ie.*Rs)./Rp)-Ie;
fd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ie.*Rs)./(a1*VT))-1))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ie.*Rs)./(a2*VT))-1))...
    -(Io3.*(Rs/a3*VT).*(exp((Vp+Ie.*Rs)./(a3*VT))-1))-(Rs/Rp)-1;
Ip=Ie-M.*(F./fd); 

FF=Iph-Io1.*[exp((Vp+Ip.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ip.*Rs)./(a2*VT))-1]...
    -Io3.*[exp((Vp+Ip.*Rs)./(a3*VT))-1]-((Vp+Ip.*Rs)./Rp)-Ip;
fdd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ip.*Rs)./(a1*VT))-1))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ip.*Rs)./(a2*VT))-1))...
    -(Io3.*(Rs/a3*VT).*(exp((Vp+Ip.*Rs)./(a3*VT))-1))-(Rs/Rp)-1;
Ip=Ip-M.*(FF./fdd); 

FFF=Iph-Io1.*[exp((Vp+Ip.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ip.*Rs)./(a2*VT))-1]...
    -Io3.*[exp((Vp+Ip.*Rs)./(a3*VT))-1]-((Vp+Ip.*Rs)./Rp)-Ip;
fddd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ip.*Rs)./(a1*VT))-1))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ip.*Rs)./(a2*VT))-1))...
    -(Io3.*(Rs/a3*VT).*(exp((Vp+Ip.*Rs)./(a3*VT))-1))-(Rs/Rp)-1;
Ip=Ip-M.*(FFF./fddd);

FFFF=Iph-Io1.*[exp((Vp+Ip.*Rs)./(a1*VT))-1]-Io2.*[exp((Vp+Ip.*Rs)./(a2*VT))-1]...
    -Io3.*[exp((Vp+Ip.*Rs)./(a3*VT))-1]-((Vp+Ip.*Rs)./Rp)-Ip;
fdddd=-(Io1.*(Rs/a1*VT).*(exp((Vp+Ip.*Rs)./(a1*VT))-1))-(Io2.*(Rs/a2*VT).*(exp((Vp+Ip.*Rs)./(a2*VT))-1))...
    -(Io3.*(Rs/a3*VT).*(exp((Vp+Ip.*Rs)./(a3*VT))-1))-(Rs/Rp)-1;
for i=1:length(Vp)
Ip(i)=Ie(i)-M*(F(i)/fd(i))-(F(i)^4/fdddd(i));
end
%%%%%%%%%%%%%%%%%//// Computing the fitness function ////%%%%%%%%%%%%%%%%%%
x_rms=sqrt((1/N)*sum((Ie-Ip).^2));
Fit=x_rms;
% 

