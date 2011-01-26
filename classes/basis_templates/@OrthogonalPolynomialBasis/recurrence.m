function[a,b] = recurrence(self, n)
% recurrence -- Evaluates the recurrence relation for the polynomial family
%
% [alpha,beta] = recurrence(n)
%
%     Computes and returns the recurrence coefficients alpha_m and beta_m for
%     each m in the array n. 

[a,b] = self.recurrence_handle(n);
