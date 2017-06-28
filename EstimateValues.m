function [pref, var, error, results] = EstimateValues(input)

max_iter = 50;
step = 1000;
results = zeros(5,max_iter-1);

x = [-400 -400 -200 -200 -100 -100 0 0 100 100 200 200 400 400]';

pref = 100;

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
        pref = -400 + 800/14*max_input;
    end
end

% Look for changes of signin simple derivative, to select points that would
% be used to define a baselin.
diff = input(2:18)-input(1:17);
first_id = 0;
last_id = 0;
for i=3:max_input-1
    if(sign(diff(i)) ~= sign(diff(i-1)))
        first_id = i;
    end
end

for i=17:-1:max_input+1
    if(sign(diff(i)) ~= sign(diff(i-1)))
        last_id = i;
    end
end

% Getting the points for the baseline
if last_id>0 && first_id>0
    lin = [2:first_id last_id:18];
    baseline =  [input(2:first_id) input(last_id:18)];
else if last_id>0
        lin = last_id:18;
        baseline =  input(last_id:18);
    else if first_id>0
            lin = 2:first_id;
            baseline =  input(2:first_id);
        end
    end
end

calculated_baseline = polyfit(lin',baseline,1);

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
    
    esti = ResponseEstimation(pref,var);
    esti = assembleFunctions(calculated_baseline,esti);
    
    error = sum(abs(input - esti'));
    results(3,iter) = error;
    
    if error > 0.1
       syms a b
       % Update pref
       old_pref = pref;
       var_c = var/(2*sqrt(2*log(2)));
       hrf = getcanonicalhrf(4,2);
       convo = convSym(exp(-1/b^2*(x-a).^2),hrf);
       convoMax = max(vpa(subs(subs(convo,a,old_pref),b,var_c)));
       convo = convo/convoMax;
       convo = convo(1:18);
       poly = polyval(calculated_baseline,1:18);
       base_in = input' - poly;
       grad_pref_eq = sum(diff(input - convo,a));
       grad_pref = vpa(subs(subs(grad_pref_eq,a,old_pref),b,var_c));
       results(4,iter) = grad_pref;
       %grad_pref = sum(-(log(2)*(pref-x)*(2.^(5-(8*(x-pref).^2)/var.^2))'*(input-2.^(-8*(x-pref).^2)/var.^2))/var.^2) ;
       pref = pref - step*grad_pref;
       fprintf('Pref : %d\n',double(pref));
       % Update var
       grad_var_eq = sum(diff(input-assembleFunctions(calculated_baseline,lin, convo'),b));
       grad_var = vpa(subs(subs(grad_var_eq,a,old_pref),b,var_c));
       results(5,iter) = grad_var;
       %grad_var = sum(-(log(2)*(old_pref-x).^2*(2.^(5-(8*(x-old_pref).^2)/var.^2))'*(input-2.^(-8*(x-old_pref).^2)/var.^2))/var.^3);
       var = var - step*grad_var;
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
subplot(3,1,1)
plot(results(1,:));
subplot(3,1,2)
plot(results(2,:));
subplot(3,1,3)
plot(results(3,:));
figure
subplot(2,1,1)
plot(results(4,:));
subplot(2,1,2)
plot(results(5,:));

end

function w = convSym(u,v)
   % CONV Nonnumeric convolution.
   % W = CONV(U,V) is the convolution of symbolic vectors U and V.
   % length(w) = length(u)+length(v)-1
   % w(k) = sum(u(j)*v(k+1-j)), j = max(1,k+1-length(v):min(k,length(u))
   
   % Form Toeplitz matrix from v.
   m = length(u);
   n = length(v);
   p = m+n-1;
   T = [v(:); zeros(m,1)]*ones(1,m);
   T = reshape(T(1:m*p),p,m);
   
   % Convolve
   w = T*u(:);
   
   % Return row vector if appropriate
   if size(u,2) > size(u,1)
      w = w.';
   end
end

function result = assembleFunctions(line_coef, conv)
    result = zeros(1,18);
    result(4:18) = polyval(line_coef,4:18);
    conv_set = 1:18;
    result(conv_set) = result(conv_set) + (1-result(conv_set)).*conv(conv_set);
    figure
    plot(0:2:34,result(1:18));
end