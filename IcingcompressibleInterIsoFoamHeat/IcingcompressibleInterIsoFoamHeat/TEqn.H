//scalar dt = runTime.deltaTValue();
Solid = alpha1*(1-f_ls);
{   
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             divUp - contErr/rho*p
          +  (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))() - contErr*K 
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
          +alpha2()/mixture.thermo2().Cv()()
        )
     == fvOptions(rho, T) - SH/IceCv
    ); 

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correctThermo();
    mixture.correct();
    
    SH = (-4.0)*rho*k_l*Solid*L*(1.0-(tanh(8.0*k_l*(T-T0)))*(tanh(8.0*k_l*(T-T0))))*fvc::ddt(T);
}

