function Fit=objfuntdnr4(Xnew,C_Iter)
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



dataset=[...
   -0.2057    0.7640
   -0.1291    0.7620
   -0.0588    0.7605
    0.0057    0.7605
    0.0646    0.7600
    0.1185    0.7590
    0.1678    0.7570
    0.2132    0.7570
    0.2545    0.7555
    0.2924    0.7540
    0.3269    0.7505
    0.3585    0.7465
    0.3873    0.7385
    0.4137    0.7280
    0.4373    0.7065
    0.4590    0.6755
    0.4784    0.6320
    0.4960    0.5730
    0.5119    0.4990
    0.5265    0.4130
    0.5398    0.3165
    0.5521    0.2120
    0.5633    0.1035
    0.5736   -0.0100
    0.5833   -0.1230
    0.5900   -0.2100
];

VT=(k*Tc)/q;
Vp=dataset(:,1);
Ie=dataset(:,2);

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

