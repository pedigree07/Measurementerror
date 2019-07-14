SEcase =  function(n1, n2, NewV, C.V, ZW, W, y1, X2, Coef.ZX, C.coef.ZX, Cum.C.coef.WS, C.Resi, Resi_zx)
{     
   X2 = as.matrix(X2)
   Ratio = n1/n2; p = dim(ZW)[2]; p2 = dim(X2)[2]
   V = (C.V+NewV)/(n1+n2)
   Cum.Resi = C.Resi + Resi_zx + t(C.coef.ZX) %*% C.V %*% C.coef.ZX + 
      t(Coef.ZX) %*% NewV %*% Coef.ZX - t(Cum.C.coef.WS) %*% (NewV + C.V) %*% Cum.C.coef.WS
   F1first = t(X2) %*% ( X2 * as.vector((X2 %*% Coef.ZX[(p-p2+1):p])^2))/(n2)
   F1second = t(X2) %*% X2 %*% Coef.ZX[(p-p2+1):p]/(n2)
   F1 = matrix(0,p,p)
   F1[(p-p2+1):p,(p-p2+1):p] = F1first - (F1second %*% t(F1second))
   Last = V * drop(Cum.Resi)/(n1+n2-p)
   first = ZW * as.vector(y1)
   second = ZW * as.vector(ZW %*%  Coef.ZX)
   zeormat = matrix(0, n1, p-p2)
   third = cbind(zeormat,W * as.vector((W %*%  Coef.ZX[(p-p2+1):p])))
   forth = cbind(zeormat, matrix(rep(F1second,n1),n1,p2,byrow=TRUE))
   Firstterm = t(first - second + third - forth) %*% (first - second + third - forth)/(n1)
   F4 = Ratio/(Ratio+1) * Firstterm + (1/(Ratio+1) * (F1*Ratio^2 + Last))
   AsymVar3 = solve(V,F4) %*% solve(V)/(n1+n2)
   return(AsymVar3 = AsymVar3)
}