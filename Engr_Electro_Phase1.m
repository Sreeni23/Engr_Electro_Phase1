%Intitialize Time
d=0.01;
mx = 100/d;
t=0:d:100;

%Get current from user
I=zeros(1,mx);
c_value=input('Enter value of current(in uA):');
c_time=input('Enter the time for which the current is applied(in ms, max:100ms):');
limit=c_time/d;
I(1:limit)=c_value;
   
%Constants
g_K = 36; g_Na = 120; g_L = 0.3; 
E_K = -12; E_Na = 115; E_L = 10.6;
Vrest = -70;
Cm = 1;

%Variables
V=zeros(1,mx);
Am=zeros(1,mx);Bm=zeros(1,mx);An=zeros(1,mx);Bn=zeros(1,mx);Ah=zeros(1,mx);Bh=zeros(1,mx);
m=zeros(1,mx);n=zeros(1,mx);h=zeros(1,mx);
gNa=zeros(1,mx+1);gK=zeros(1,mx+1);

%Initial Conditions 
V0 = 0;
Am0 = 0.1*((25-V0)/(exp((25-V0)/10)-1));
Bm0 = 4*exp(-V0/18);
An0 = 0.01*((10-V0)/(exp((10-V0)/10)-1));
Bn0 = 0.125*exp(-V0/80);
Ah0 = 0.07*exp(-V0/20);
Bh0 = 1/(exp((30-V0)/10)+1);

m(1) = Am0/(Am0+Bm0);
n(1) = An0/(An0+Bn0);
h(1) = Ah0/(Ah0+Bh0);

%Calculations
    for i=1:mx
        Am(i) = 0.1*((25-V(i))/(exp((25-V(i))/10)-1));
        Bm(i) = 4*(exp(-V(i)/18));
        An(i) = 0.01*((10-V(i))/(exp((10-V(i))/10)-1));
        Bn(i) = 0.125*exp(-V(i)/(80));
        Ah(i) = 0.07*exp(-V(i)/(20));
        Bh(i) = 1/(exp((30-V(i))/10)+1);
        
        gNa(i) = (m(i)^3)*g_Na*h(i);
        gK(i) = (n(i)^4)*g_K;
        I_Na = gNa(i)*(V(i)-E_Na);
        I_K = gK(i)*(V(i)-E_K);
        I_L = g_L*(V(i)-E_L);
             
        I_ion = I(i)-I_K-I_Na-I_L;   
        
        V(i+1) = V(i)+d*(I_ion/Cm);
        m(i+1) = m(i)+d*(Am(i)*(1-m(i))-Bm(i)*m(i));
        n(i+1) = n(i)+d*(An(i)*(1-n(i))-Bn(i)*n(i));
        h(i+1) = h(i)+d*(Ah(i)*(1-h(i))-Bh(i)*h(i));
   end
gNa(mx+1)=gNa(mx);
gK(mx+1)=gK(mx);

Vm=V+Vrest;
figure(1);
subplot(1,2,1)
plot(t,Vm)
title('V_membrane vs Time')
ylabel('Membrane Voltage(mV)')
xlabel('Time(ms)')
subplot(1,2,2)
plot(t,gNa,t,gK)
title('Conductance of K and Na vs Time')
ylabel('Conductance(mS/cm2)')
xlabel('Time(ms)')
legend('gNa','gK')
clear;  