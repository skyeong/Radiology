function [fa fa2] = CalcProperties(X)

if size(X,2)~=3,
    X = X';
end

[coeff, score, latent, tsquared, explained] = pca(X,'NumComponents',3);
% [U,S,V] = svd(X);
% lambda = [S(1,1) S(2,2) S(3,3)];
% lambda = sqrt(lambda);
lambda = latent;
muL = mean(lambda);
fa = sqrt(3/2) * sqrt(sum((lambda-muL).^2)) / sqrt(sum(lambda.*lambda));


L2 = var(X);
muL2 = mean(L2);
fa2 = sqrt(3/2) * sqrt(sum((L2-muL2).^2)) / sqrt(sum(L2.*L2));


% SS = var(X);
% sortedSS = sort(SS,'descend');
% vmm = sortedSS(1)./sortedSS(2);
