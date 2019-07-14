Datageneration = function(n1, n2, mu,umu, Sigma, SigmaU, model.error, q, theta)
{
   # block 1 data:
   dmat1 = mvrnorm(n = n1, mu = mu, Sigma)
   W_1 = dmat1[,-c(1:q)] + mvrnorm(n = n1, rep(umu,dim(SigmaU)[2]), SigmaU)
   Z_1 = dmat1[,1:q]
   Z_1 = cbind(1, Z_1)
   X_1 = dmat1[,-c(1:q)]
   y1 = cbind(1, dmat1) %*% theta + rnorm(n = n1,0,model.error)
   
   # block 2 data:
   dmat2 = mvrnorm(n = n2, mu= mu, Sigma)
   X_2 = dmat2[,-c(1:q)]
   Z_2 = dmat2[,1:q]
   Z_2 = cbind(1, Z_2)
   y2 = cbind(1, dmat2) %*% theta + rnorm(n = n2,0,model.error)
   list(y1 = y1, y2 = y2, W_1 = W_1, Z_1 = Z_1, X_2 = X_2, Z_2 = Z_2, X_1 = X_1)
}