% run_sumulation_v4_1.m, V. Ziemann, 231212
% uses the forward signals of the currents
% x=[Vreal;Vimag];   
% u=[Ireal_f,Iimag_f];
clear all; % close all
Nforget=1000000;    % forgetting horizon
alpha=1-1/Nforget;
Npulse=1000;    % negative values turn off pulses
R=1;            % shunt impedance, ensure current and voltage is normalized
sigp=0.0001;    % process noise level
sigm=0.001;     % measurement noise level
omega0=2*pi*1e9;
QE=1e6;                   % external-Q
%QE=8e8;
Q0=1e9;                   % unloaded Q0
QL=1/(1/QE+1/Q0);         % loaded-Q
omega12=omega0/(2*QL);    % bandwidth
omegaE=omega0/QE;         
domega=omega12/2;        % set detuning to half the bandwidth
%domega=0;
if QE<2e8    
  dt=1e-4;        % sample time at 10kHz
  Niter=60000;     % number of iterations
else
  dt=1e-3;
  Niter=10000;     % number of iterations
end
q0=[omega12*dt,domega*dt];   % bandwidth and detuning
qfit=[omegaE*dt,domega*dt,omega0*dt/Q0];
bandwidth=q0(1);
F0=[-q0(1),-q0(2);q0(2),-q0(1)];  % eq. 4, F0=Abar 
Areal=eye(2)+F0;            % eq. 6
Breal=R*omega12*dt*eye(2);  % eq. 1, for generator current, used to drive cavity
if 0 %..............................calibration, set to "1" for Fig.3
  tic
  disp("Forcing detuning to zero!"); beep;
  Areal(1,2)=0; Areal(2,1)=0; % make sure that NO detuning == on-resonance
  alpha=1;
  Niter=1000000;
  xp=sigp*randn(2,1);         % initialize cavity voltage reading, process noise
  x=sigm*randn(2,1);          % measured voltage inside cavity
  u=[1;0];
  for iter=1:Niter % iterate to reach equilibrium
    xpnew=Areal*xp+Breal*u+sigp*randn(size(x)); % cavity dynamics
    xnew=xpnew+sigm*randn(size(x)); 
    xp=xpnew;       % update process voltage in the cavity
  end
  qhat=0;  
  data=zeros(Niter,2);
  pp=1;  % PT for one degree of freedom
  for iter=1:Niter % sysid to determine R
    xpnew=Areal*xp+Breal*u+sigp*randn(size(x)); % cavity dynamics
    xnew=xpnew+sigm*randn(size(x)); 
    u=[1;0]; %+sigm*randn(2,1);
    tmp2=1/(alpha+pp*u'*u);    % sysid for R
    ppnew=tmp2*pp;                        
    qhat=tmp2*(alpha*qhat+pp*u'*xnew);
    xp=xpnew;       % update process voltage in the cavity
    pp=ppnew;
    data(iter,1)=qhat;
    data(iter,2)=pp;
  end
  mm=1:Niter;         % xaxis for plots
  figure(1); clf; %subplot(2,1,1);
    loglog(mm,abs(data(:,1)-1),'k','LineWidth',2);
  %loglog(mm,abs(data(:,1)-1),'k',mm,sigm*sqrt(data(:,2)),'r--','LineWidth',2);
  xlim([1,Niter]); ylim([6e-7,1])
  xlabel('Iterations'); ylabel('\deltaR');
  set(gca,'FontSize',16);
  figure(2); %subplot(2,1,2);
  loglog(mm,data(:,2),'k','LineWidth',2);
  xlim([1,Niter]);
  xlabel('Iterations'); ylabel('P_T'); 
  set(gca,'FontSize',16);
  R=qhat;
  toc
  return
end  % end of calibration
P=eye(3);                   % initial value of pT
qhat=zeros(3,1);            % initial parameter estimate
xp=sigp*randn(2,1);         % initialize cavity voltage reading, process noise
x=sigm*randn(2,1);          % measured voltage inside cavity
data=zeros(Niter,10);       % storage for later plotting
uset=[0;0];                 % generator is off at start
if Npulse>0, disp(['Pulsing with period ',num2str(Npulse)]); beep;  end
eee=0e-2;    % error on R and voltage scale (anti-correlated)
scal=1;      % scale factor that affects voltage and R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for iter=1:Niter                 % main iteration loop
  %if iter==100, uset=[1;0]; end  % first pulse
  if Npulse>0                    % negative Npulse turns off later pulses
    if mod(iter,Npulse)==1       % pulsed operation
      if uset(1)==0, uset=[1;0]; else, uset=[0;0]; end
    end
  end
  %.........................................cavity dynamics
  u=uset;
  xpnew=Areal*xp+Breal*u+sigp*randn(size(x)); % eq. 2, xp=V
  xnew=scal*xpnew+sigm*randn(size(x));        % xnew=V', add measurement noise
  %...................................system identification
  %up=u*omega12/omegaE;             % forward current I^+ 
  up=u*QE/(2*QL);                   % up=forward current I^+ 
  y=xnew-x;                                       
  G=[-0.5*x(1)+up(1)*R*scal,-x(2),-0.5*x(1);-0.5*x(2)+up(2)*R*scal,x(1),-0.5*x(2)];
  tmp2=eye(3)-P*G'*inv(alpha*eye(2)+G*P*G')*G;    
  Pnew=tmp2*P/alpha;                              
  qhat=tmp2*(qhat+P*G'*y/alpha);
  %.........................................................
  xp=xpnew;       % update process voltage in the cavity
  x=xnew;         % remember measured voltage
  P=Pnew;
  %................................save for later plotting
  data(iter,1)=x(1);       % normalized voltages, measured
  data(iter,2)=x(2);
  data(iter,3)=u(1);       % generator currents
  data(iter,4)=u(2);
  data(iter,5)=qhat(1);    % omegaE*dt
  data(iter,6)=qhat(2);    % detuning
  data(iter,7)=qhat(3);    % omega00*dt
  data(iter,8)=P(1,1);     % covariance matrix diagonals
  data(iter,9)=P(2,2);
  data(iter,10)=P(3,3);  
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
  if QE<2e8   % over-coupled
    semilogy(mm,fac*data(:,5),'k',mm,fac*data(:,6),'r',mm,abs(fac*data(:,7)),'b--','LineWidth',2);
    ylim([0.01,1500])
  else
    plot(mm,fac*data(:,5),'k',mm,fac*data(:,6),'r',mm,fac*data(:,7),'b--','LineWidth',2);
  end
  xlabel('Iterations'); ylabel('f_{E},\Deltaf, f_{Q0}   [Hz]');
  legend('f_{E}','\Deltaf','f_{Q0}');
  set(gca,'FontSize',16);
  %fit_at_end=[data(end,5),data(end,6),data(end,7)]*fac
end

if 1   % evolution of the Q-values
  figure(3); clf
  QQE=omega0*dt./data(:,5);
  QQ0=omega0*dt./data(:,7);
  if QE< 2e8
    semilogy(mm,QQE,'k',mm,QQ0,'b--','LineWidth',2)
    ylim([4e5,8e9])
  else
    plot(mm,QQE,'k',mm,QQ0,'b--','LineWidth',2)
    ylim([6e8,1.2e9])
  end
  xlabel('Iterations'); ylabel('Q_E, Q_0');
  legend('Q_E','Q_0');
  set(gca,'FontSize',16);
  %Q_at_end=omega0*dt./[data(end,5),data(end,7)]*1e-6
  disp(['QE = ',num2str(1e-6*omega0*dt/data(end,5)),'E6,  Q0 = ',num2str(1e-6*omega0*dt/data(end,7)),'E6'])
end
% if QE<2e8
%   figure(5); clf
%   QQ0=abs(omega0*dt./data(:,7));
%   semilogy(mm,QQ0,'b--','LineWidth',2)
%   xlabel('Iterations'); ylabel('Q_0');
%   set(gca,'FontSize',16);
% end

if 1   % error bars
  figure(4); clf %   plot PT
  eb1=sqrt(data(:,8))*sigm;
  eb2=sqrt(data(:,9))*sigm;
  eb3=sqrt(data(:,10))*sigm;
  % semilogy(mm,abs(eb1./data(:,5)),'k',mm,abs(eb2./data(:,6)),'r:', ...
  %  mm,abs(eb3./data(:,7)),'b','LineWidth',2);
  % legend('Q_E','\Delta\omega','Q_0');
  semilogy(mm,abs(eb1./data(:,5)),'k',mm,abs(eb3./data(:,7)),'b--','LineWidth',2);
  ylim([5e-5,2]);
  xlim([1,Niter]); 
  xlabel('Iterations'); ylabel('rel. uncertainty')
  legend('Q_E','Q_0');
  set(gca,'FontSize',16);
end

if 0   % estimation error
  figure(5); clf                  % Plot estimation error
  eb1=sqrt(sigm^2+2*sigp^2)*sqrt(data(:,7)); % empirical error bar of q(1)
  %eb1=sigm*sqrt(pt(mm));
  data(:,1)=abs(data(:,5)-q0(1)); % estimation error of q(1)
  loglog(mm,data(:,1),'k',mm,eb1,'b-.','LineWidth',2);
  xlim([1,Niter]); ylim([7e-7,2*max(data(:,5))])
  xlabel('Iterations'); ylabel('Estimation error |a_T(1)|')
  legend('Simulation','Errorbars');
  set(gca,'FontSize',16);
end
