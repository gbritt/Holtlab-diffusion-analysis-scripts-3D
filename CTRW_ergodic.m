function phi = CTRW_ergodic(a,eps)


phi = ( (gamma(1+a)^(1/a)) ./ (a * eps.^(1+1/a)) ) .* sqrt(a/(2*pi)) .* eps.^(3/(2*a)) / gamma(1+a)^(3/(2*a)).* ...
    exp(- (a/2) * eps.^(1/a) / gamma(1+a)^(1/a) );



end