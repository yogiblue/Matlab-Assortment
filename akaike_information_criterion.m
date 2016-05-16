for p=1:20
    [A, E(p)]=arburg(x,p);
    AIC(p) = log(E(p)) + 2*p/length(x);
end
plot(AIC)
