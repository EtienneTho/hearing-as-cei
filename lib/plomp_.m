function plomp = plomp_(f1,f2)
    fmin=min(f1,f2);
    fmax=max(f1,f2);
    s=0.24/(0.0207*fmin+18.96);
    plomp = exp(-3.5*s*(fmax-fmin))-exp(-5.75*s*(fmax-fmin));
end