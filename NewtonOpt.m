%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project 4: Newton Optimization
% MA332 - 01: Introduction to Comp. Sci.
% Dr. Leader
% 30 OCT 2020
%
% Using Newtons method to solve unconstrained optimization problems
% utilizing the first order(gradient) condition of optimality and the
% second order(Hessian) condition of optimality. In addition, we define the
% search direction (minimizing -f'(x)) and the step size (Lagrange
% polynomials alpha). Using Newton Optimization, we can confidently
% computationally model real world problems/functions.
%
% Author: Raymond Becerra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbest, fbest, itrcnt, stat] = NewtonOpt(f, Df, Df2, x0, gradTol, xTol, itrBound, s)

stat    = 1;
itrcnt  = 1;
xbest   = x0;
stop    = false;

while stop == false
    
    grad    = ( (feval(Df2, xbest)) \ (-feval(Df, xbest)) );
    
    if s == 1
        
        alpha   = 1;
    
    elseif s == 2
        
        % Assigning Lagrange polynomial values
        a1      = 0;
        a2      = 1;
        a3      = 2;
        
        % Assigning taylor polynomial functions
        F1    = feval(f, xbest);
        F2    = feval(f, xbest + grad);
        F3    = feval(f, xbest + 2*grad);
        
        if abs( (2*(a2 - a3))*F1 + (2*(a3 - a1))*F2 + (2*(a1 - a2))*F3 ) <= (10^-10);
            
            alpha   = 1;
            
        else
            % Using equation 4.7 from the textbook to calculate value of
            % alpha for lagrange polynomial a = 2 and a valid denominator
            % tolerance. This alpha is used to comupte the step size for
            % the search direction
            num     = ( (a2^2) - (a3^2) )*F1 + ( (a3^2) - (a1^2) )*F2 + ( (a1^2) - (a2^2) )*F3;
            den     = ( (2*(a2 - a3))*F1 + (2*(a3 - a1))*F2 + (2*(a1 - a2))*F3 );
            alpha   = num / den;
            
        end
        
    elseif s == 3
            
        % Using the bisection method from Chap.1 to compute an alpha value
        % for the computation of the step size for search direction.
        bisectAlpha     = @(aHat) (feval(Df, xbest + aHat*grad)' * (grad) );
        n               = abs( round( ( -8 * log(10) ) / log(2) ) );
        aL_limit        = 0;
        aR_limit        = 1;
        aLeft           = feval(bisectAlpha, aL_limit);
        aRight          = feval(bisectAlpha, aR_limit);
        
    for run = 1:n
        
        aMid    = (aL_limit + aR_limit) / 2;
        fMid    = feval(bisectAlpha, aMid);
        
    if sign(fMid) == 0
        
        aR_limit    = aMid;
        
        %break;
        
    elseif ( sign(fMid) == sign(aLeft) )
        
        aL_limit    = aMid;
        aLeft       = fMid;
        
    elseif ( sign(fMid) == sign(aRight) )
        
        aR_limit    = aMid;
        aRight      = fMid;
        
    end
    
    end
    
        alpha = aR_limit;
        
    elseif s == 4
        % Using equation 4.6 to calculate alpha for a second order taylor
        % polynomial using gradient and the Hessian.
        alpha   = -( (feval(Df, xbest)' * grad) / (grad' * feval(Df2, xbest)' * grad) );
        
    else
        
        alpha   = 1;
        
    end
    
    xNew    = xbest + alpha*grad;
    Df_norm = norm( feval(Df, xNew) );
    
    if Df_norm <= gradTol
        
        stop    = true;
        stat    = 0;
        
    end
    
    if norm(xNew - xbest) <= xTol
        
        stop    = true;
        
    end
    
    if itrcnt >= itrBound
        
        stop    = true;
        
    end
    
    xbest   = xNew;
    itrcnt = itrcnt + 1;
end

fbest   = feval(f,xbest);

end
