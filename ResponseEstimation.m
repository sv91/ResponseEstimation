function [result] = ResponseEstimation(pref, var,baseline)
warning('off','all');
%% Define curve
toPlot = false;
% Parameter
x_pref = pref;
FWHM = var;

var = FWHM/(2*sqrt(2*log(2)));
y = zeros(801,1);
x= -400:400;

for i=1:801
    y(i) = baseline + (1-baseline) * exp(-1/var^2*(x(i)-x_pref)^2);
end

if(toPlot)
    figure(1)
    subplot(2,2,1);
    plot(x,y);
end

%% Discrete values

x_disc = -400:50:400;
s_disc = size(x_disc,2);
y_disc = zeros(801,1);

for i=1:s_disc
    y_disc(x_disc(i)+401) = baseline + (1-baseline)*exp(-1/var^2*(x_disc(i)-x_pref)^2);
end

if(toPlot)
    subplot(2,2,2);
    plot(x,y_disc);
end

%% Experience
t_disc = 0:4:35;
x_disc = [-400 -200 -100 0 100 200 400 NaN NaN];
y_disc = baseline + (1-baseline)*exp(-1/var^2*(x_disc-x_pref).^2);

max_y = max(y_disc);

if(toPlot)
    subplot(2,2,3);
    plot(t_disc,y_disc);
end
%% Transform to a 36 points matrix
t_f = 0:2:35;
y_f = zeros(36,1);
for i=1:36
  if(mod(i,4)==1)
      y_f(i) = y_disc((i-1)/4+1);
  end
  if(mod(i,4)==3)
      y_f(i) = y_disc((i-3)/4+1);
  end
   y_f(i) = y_disc(floor((i-1)/4)+1);
end

%% Convolve
% resample to desired TR
y_f = interp1((0:length(y_f)-1),y_f,0:2:(length(y_f)-1),'pchip');
hrf = getcanonicalhrf(4,2);
convo = conv(y_f,hrf);
convo = convo/max(convo);

% [~,max_val] = max(convo);
% for i=max_val:18
%    convo(i)= baseline+(1-baseline)*convo(i); 
% end
if(toPlot)
    subplot(2,2,4);
    plot(t_f,convo(1:18));
end

warning('on','all');
result = convo(1:18);