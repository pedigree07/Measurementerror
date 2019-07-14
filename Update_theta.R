OnlineMeasearment = function(y1, y2, Z_1, W_1, Z_2, X_2, X_1)
{
   ZW = cbind(Z_1, W_1); ZX = cbind(Z_2, X_2); X1Z = cbind(Z_1, X_1); Full = rbind(X1Z,ZX);
   y = c(y1,y2); n_1 = dim(ZW)[1]; n_2 = dim(ZX)[1]; p = dim(ZW)[2]
   ######Correct biased Estimates
   V =  t(ZW) %*% ZW
   NewV =  t(ZX) %*% ZX
   Coef.ZW = solve( V ) %*% (t(ZW) %*% y1)
   Coef.ZX = solve( NewV ) %*% (t(ZX) %*% y2)
   Resi_zw = sum( (y1 - ZW %*% Coef.ZW)^2 ) 
   Resi_zx = sum( (y2 - ZX %*% Coef.ZX)^2 ) 
   var_zx = solve( NewV ) * Resi_zx/( n_2 - p )
   W1W1 = t(W_1)%*%W_1; Z1Z1 = t(Z_1)%*%Z_1; X2X2 = t(X_2)%*%X_2* n_1/n_2
   Z1W1 = t(Z_1)%*%W_1; W1Z1 = t(Z1W1); SigU = W1W1 - X2X2;
   C.V = rbind( cbind(Z1Z1, Z1W1), cbind(W1Z1, X2X2) )
   C.coef.ZX = solve(C.V) %*% (t(ZW) %*% y1) 
   ######Update Estimates
   Delta =  Coef.ZW - C.coef.ZX
   Umat = matrix(0, nrow = p, ncol = p)
   Umat[ (dim(Z_1)[2] + 1):p, (dim(Z_1)[2] + 1):p] = W1W1 - X2X2
   First_T = t(Delta) %*% C.V %*% Delta
   C.Resi = Resi_zw - First_T - t(Coef.ZW) %*%  Umat %*% Coef.ZW
   Cum.C.coef.WS = solve( C.V + (NewV) ) %*%
      ( C.V %*% C.coef.ZX + (NewV) %*% Coef.ZX)
   list(Coef.ZX = drop(Coef.ZX), C.coef.ZX = drop(C.coef.ZX), Cum.C.coef.WS = drop(Cum.C.coef.WS),
        C.Resi = C.Resi, Resi_zx = Resi_zx, V = V, C.V = C.V, NewV = NewV, ZW = ZW)
}