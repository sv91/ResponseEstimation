function [pref, var, error, results] = EstimateValues(input, varargin)


%% Managing the inputs
if nargin > 0
   for i=1:2:nargin
       switch varargin{i}
           case 'baseline'
               fit_bl        = varargin{i+1};
           case 'maxiter'
               max_iter      = varargin{i+1};
           case 'last'
               copy_last     = varargin{i+1};
           case 'error'
               error_limit   = varargin{i+1};
           case 'step'
               step          = varargin{i+1};
           case 'baseline_step'
               b_step        = varargin{i+1};
           otherwise
               fprintf('Unknown argument "%s"\n',varargin{i});
       end
   end
end

%% Setting to default all the unset variables.
if notDefined('fit_bl')   
    fit_bl      = false;
end

if notDefined('max_iter')   
    max_iter    = 50;
end

if notDefined('copy_last')   
    copy_last   = false;
end

if notDefined('error_limit')   
    error_limit = 0.1;
end

if notDefined('step')   
    step = 400;
end

%%
results = zeros(7,max_iter-1);

x = [-400 -400 -200 -200 -100 -100 0 0 100 100 200 200 400 400]';
syms a b c d;      
% Create the equations
hrf = getcanonicalhrf(4,2);
convo = convSym(c + (1-c)*exp(-1/b^2*(x-a).^2),hrf);
% grad_base_eq = sum(diff(abs(input - convo(1:18)/d),c));
grad_pref_eq = sum(diff((input - convo(1:18)/d).^2,a));
grad_var_eq = sum(diff((input - convo(1:18)/d).^2,b));


pref = 100;
baseline = 0;

% Find a default value for the preffered SOA, as the equivalent SOA to the
% index of the max value of the input.
[~,max_input] = max(input);
if max_input < 3
   if max_input >= 2
       pref = -400;
   else 
       return;
   end
else
    if max_input > 17
        return;
    else
        pref = -400 + 800/14*(max_input-3);
    end
end
baseline = mean(input(setdiff(3:18,max_input-1:max_input+1)));
clear max_input

var = 100;
error = 1000;
iter = 1;
esti = zeros(size(input));

while iter < max_iter && error > 0.1
    if(mod(iter,10)==0)
       fprintf('Iteration %d\n',iter);
    end
    results(1,iter) = pref;
    results(2,iter) = var;
    results(3,iter) = baseline;
    
    esti = ResponseEstimation(pref,var,baseline);
    %esti = assembleFunctions(calculated_baseline,esti);
    plot_input = [input; input(1)];
    plot_esti = [esti esti(1)];
    figure
    plot(0:2:36,plot_esti,0:2:36,plot_input);
    clear plot_input plot_esti
    
    error = sum((input - esti').^2)
    results(4,iter) = error;
    
    if error > 0.1
       % Update pref
       var_c = var/(2*sqrt(2*log(2)));
       old_baseline = baseline;
       old_pref = pref;

       convoMax = max(vpa(subs(subs(subs(convo,a,pref),c,baseline),b,var_c)));
       
       % Updating the baseline by using grid search and binomial search
       range_base = (1/10)^(1+(iter-1)/max_iter);
       temp_baseline = baseline;
       for i=1:3
           range_base = range_base/2;
           base_p =  temp_baseline + range_base;
           base_m =  temp_baseline - range_base;
           esti_p = ResponseEstimation(pref,var,base_p);
           error_p= sum((input - esti_p').^2);
           esti_m = ResponseEstimation(pref,var,base_m);
           error_m= sum((input - esti_m').^2);
           if error_p<error_m
               temp_baseline = base_p;
               error_temp = error_p;
           else
               temp_baseline = base_m;
               error_temp = error_m;
           end
           clear base_p base_m esti_p error_p esti_m error_m
       end
       if error_temp < error
           baseline = temp_baseline;
       end
       clear range_base temp_baseline
       
       
%        grad_base = vpa(subs(subs(subs(subs(grad_base_eq,a,pref),d,convoMax),c,old_baseline),b,var_c));
%        results(5,iter) = grad_base;
%        baseline = baseline - step_b*grad_base;
%        clear grad_base
%        if baseline>1
%            baseline = 1;
%        end
%        if baseline<0
%            baseline = 0;
%        end
       fprintf('Baseline : %d\n',double(baseline));
       
       grad_pref = vpa(subs(subs(subs(subs(grad_pref_eq,a,pref),d,convoMax),c,old_baseline),b,var_c));
       results(6,iter) = grad_pref;
       pref = pref - step*grad_pref;
       clear grad_pref
       if pref>500
           pref = 500;
       end
       if pref<-500
           pref = -500;
       end
       fprintf('Pref : %d\n',double(pref));
       % Update var
       grad_var = vpa(subs(subs(subs(subs(grad_var_eq,a,old_pref),d,convoMax),c,old_baseline),b,var_c));
       results(7,iter) = grad_var;
       var = var - step*grad_var;
       if var<1
           var = 1;
       end
       clear grad_var
       fprintf('var : %d\n',double(var));
    end
    iter = iter +1;
end

figure
subplot(2,1,1)
plot(0:2:35,input);
subplot(2,1,2)
plot(0:2:35,esti);
fprintf('Pref: %d     Var: %d ',double(pref),double(var));
figure
subplot(4,1,1)
plot(results(1,:));
subplot(4,1,2)
plot(results(2,:));
subplot(4,1,3)
plot(results(3,:));
subplot(4,1,4)
plot(results(4,:));
figure
subplot(3,1,1)
plot(results(5,:));
subplot(3,1,2)
plot(results(6,:));
subplot(3,1,3)
plot(results(7,:));

end