function [kl_dist]=KL_dirichlet(alpha_1,alpha_2)

t1= gammaln(sum(alpha_1))-gammaln(sum(alpha_2));

t2= sum(gammaln(alpha_2)-gammaln(alpha_1));

t3= (alpha_1-alpha_2)'*(psi(alpha_1)-psi(sum(alpha_1)));

kl_dist=t1+t2+t3;

end