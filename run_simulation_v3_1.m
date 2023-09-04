% run_sumulation_v3_1.m, V. Ziemann, 230902
% uses the forward signals of the currents
% x=[Vreal;Vimag];   
% u=[Ireal_f,Iimag_f];
clear all; % close all
Niter=10000;     % number of iterations
Nforget=100;    % forgetting horizon
alpha=1-1/Nforget;
Npulse=-40000;  % negative values turn off pulses
R=1;            % shunt impedance, ensure current and voltage is normalized
sigp=0.0001;    % process noise level
sigm=0.001;     % measurement noise level
dt=1e-7;        % sample time at 10 MHz
%omega12=2e5;    %=1/tau=filling time=1/10 microseconds
%domega=-1e5;    %=cavity detuning =10kHz
if 1    % Spoke parameters
  omega0=2*pi*352e6;        % resonant frequency
  QE=1.8e5;                 % external-Q
  QL=QE;                    % loaded-Q, if sc same as QE
  omega12=omega0/(2*QL);    % bandwidth
  omegaE=omega0/QE;         
  domega=-omega12/2;         % set detuning to half the bandwidth
end
q0=[omega12*dt,domega*dt]   % bandwidth and detuning
bandwidth=q0(1);
F0=[-q0(1),-q0(2);q0(2),-q0(1)];  % eq. 4, F0=Abar 
Areal=eye(2)+F0;            % eq. 6
Breal=R*omega12*dt*eye(2);  % eq. 1, for generator current, used to drive cavity
BrealE=R*omegaE*dt*eye(2);  % eq. 3, for forward current used in sysid
pp=1;                       % initial value of pT
qhat=zeros(2,1);            % initial parameter estimate
xp=sigp*randn(2,1);         % initialize cavity voltage reading, process noise
x=sigm*randn(2,1);          % measured voltage inside cavity
data=zeros(Niter,7);        % storage for later plotting
uset=[0;0];                 % generator is off at start
if Npulse>0, disp(['Pulsing with period ',num2str(Npulse)]); beep;  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for iter=1:Niter                 % main iteration loop
  if iter==100, uset=[1;0]; end  % first pulse
  if Npulse>0                    % negative Npulse turns off later pulses
    if mod(iter,Npulse)==0       % pulsed operation
      if uset(1)==0, uset=[1;0]; else, uset=[0;0]; end
    end
  end
  if 1   % bandwidth change, Figure 5
    factor=2;               % magnitude of step in bandwidth
    if iter==Niter/2        % start half-way of Niter
      q0(1)=factor*q0(1);   % increase bandwidth by factor
      omega12=omega12*factor;
      Areal=eye(2)+[-q0(1),-q0(2);q0(2),-q0(1)];
      Breal=Breal*factor;
      QL=QL/factor;   
    end
    if 0 %iter==Niter/2     % TURNED OFF, stop at half of Niter
      q0(1)=q0(1)/factor;   % reduce bandwidth by factor
      omega12=omega12/factor;
      Areal=eye(2)+[-q0(1),-q0(2);q0(2),-q0(1)];
      Breal=Breal/factor;
      QL=QL*factor;
    end
  end
  if 0 % periodic perturbation (1=1 kHz, 20=20 kHz), Figure 3 and 4
    q0(2)=0.5*bandwidth*sin(20*6.2832e-4*iter);  
    Areal=eye(2)+[-q0(1),-q0(2);q0(2),-q0(1)];
  end
  %.........................................cavity dynamics
  u=uset;
  xpnew=Areal*xp+Breal*u+sigp*randn(size(x)); % eq. 4, xp=V
  xnew=xpnew+sigm*randn(size(x));             % eq. 5, x=V', add measurement noise
  %...................................system identification
  %up=u*omega12/omegaE;             % forward current I^+ 
  up=u*QE/(2*QL);                   % eq.2
  y=xnew-x-BrealE*up;               % eq.11, right
  vv2=x'*x; 
  tmp=alpha/(alpha+pp*vv2);         % bracket in eq.20
  qhat=tmp*(qhat+[-x(1)*y(1)-x(2)*y(2);-x(2)*y(1)+x(1)*y(2)]*pp/alpha); % eq.21
  pp=tmp*pp/alpha;                  % eq. 20
  xp=xpnew;       % update process voltage in the cavity
  x=xnew;         % remember measured voltage
  %................................save for later plotting
  data(iter,1)=x(1);       % normalized voltages, measured
  data(iter,2)=x(2);
  data(iter,3)=u(1);       % generator currents
  data(iter,4)=u(2);
  data(iter,5)=qhat(1);    % bandwidth
  data(iter,6)=qhat(2);    % detuning
  data(iter,7)=pp;         % p_T
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
mm=1:Niter;         % xaxis for plots
mm2=Niter/2:Niter;  % just the second half of the data

if 1    % plot of voltages and currents
  figure(1)
  subplot(2,1,1);
  plot(mm,data(:,1),'k',mm,data(:,2),'r','LineWidth',2);
  xlabel('Iterations'); ylabel('v_r, v_i');
  legend('v_r','v_i')
  %ylim([-4,10.2])
  set(gca,'FontSize',16);
  subplot(2,1,2);
  plot(mm,data(:,3),'k',mm,data(:,4),'r','LineWidth',2);
  xlabel('Iterations'); ylabel('i_r, i_i');
  legend('i_r','i_i')
  set(gca,'FontSize',16);
end

if 1   % evolution of fit parameters
  figure(2); clf
  fac=1/(2*pi*dt);   % convert q-units to Hz
  plot(mm,fac*data(:,5),'k',mm,fac*data(:,6),'r','LineWidth',2);
  xlabel('Iterations'); ylabel('f_{12},\Deltaf     [Hz]');
  legend('f_{12}','\Deltaf');
  ylim([-700,2300])
  %xlim([4950,5350])
  set(gca,'FontSize',16);
  rms_second_half=std(data(mm2,[5,6]))
end

if 1   % PT
  figure(3); clf %   plot PT
  semilogy(mm,data(:,7),'k','LineWidth',2);
  xlim([1,Niter]);  
  xlabel('Iterations'); ylabel('p_T')
  set(gca,'FontSize',16);
  asymp=1/(Nforget*vv2);
  hold on; semilogy([1,Niter],[asymp,asymp],'r--');
  error_bar_q=sqrt(asymp)*sqrt(sigm^2+2*sigp^2)
end

if 0   % estimation error
  figure(4); clf                  % Plot estimation error
  eb1=sqrt(sigm^2+2*sigp^2)*sqrt(data(:,7)); % empirical error bar of q(1)
  %eb1=sigm*sqrt(pt(mm));
  data(:,1)=abs(data(:,5)-q0(1)); % estimation error of q(1)
  loglog(mm,data(:,1),'k',mm,eb1,'b-.','LineWidth',2);
  xlim([1,Niter]); ylim([7e-7,2*max(data(:,5))])
  xlabel('Iterations'); ylabel('Estimation error |a_T(1)|')
  legend('Simulation','Errorbars');
  set(gca,'FontSize',16);
end
