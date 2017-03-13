function max_error = jacobi(Niter)

%------------------------------------------ %
% Based on the provided code for the jacobi %
% method on the unit square                 %
%-------------------------------------------%
  M = 30; 
  N = 30;
  
  len = 5; %TODO: find best value
  
  dr = (len-1) / M;
  r_vals = dr * (0:M) + 1;
  dth = 2 * pi/ N;
  theta_vals = dth * (0:N - 1);
  
  uold = zeros(M, N);
  uold(1, :) = (1-len)^3 * sin(theta_vals);

  unew=uold;
  
  % Source
  
  s = zeros(M, N);
  for i=1:M-1
      for j=1:N
          sinVal = sin(theta_vals(j));
          r = r_vals(i);
          diff = r - len;
          s(i, j) = sinVal*(6*diff ...
                           + 3*diff*diff/r ...
                           - diff^3/(r*r));
      end
  end
  
  actual = zeros(M, N);
  
  for i=1:M-1
      for j=1:N 
          r = r_vals(i);
          theta = theta_vals(j);
          actual(i, j) = (r-len)^3*sin(theta);
      end
  end
  
  % iteration
  for k=1:Niter
    for i=2:M-1
        r = r_vals(i);
            
        a = 1/(dr*dr);
        b = 1/(2*r*dr);
        
        m1 = a + b;
        m2 = a - b;
        n1 = 1/((r*r)*(dth*dth));
        u_coeff = 2*(a + n1); % Not by definition, just avoiding recomputing
        
        
        for j=1:N
            j_inc = mod(j, N) + 1; % Circle edge cases
            j_dec = mod(j-2, N) + 1; 
            
            unew(i,j) = (m1*uold(i+1, j) ...
                         + m2*uold(i-1, j) ...
                         + n1*(uold(i, j_inc) + uold(i, j_dec)) ...
                         - s(i, j)) ...
                / u_coeff;
            
        end
    end
  % For Jacobi, find all the new iterates before updating u
  % Store two arrays
    uold = unew;
  end

  u = unew;
  errors = actual - u
  max_error = mean(mean(abs(actual - u)));