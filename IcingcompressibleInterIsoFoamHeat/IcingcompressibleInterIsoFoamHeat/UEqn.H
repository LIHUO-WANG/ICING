k_l = 8.0*Pi*0.0478*0.0478*0.0478*T0*T0*(2.0*T-Tset)*f0*exp((-4.7e-7)*mag(E)/nE)/(3.0*2.87*2.87*Kb*Tset*Tset*Tset*T*T*T1);
f_ls = 0.5 + 0.5*tanh(8.0*k_l*(T-T0));
Sm = alpha1*A_mush*(1-f_ls)*(1-f_ls)/(f_ls*f_ls*f_ls+little)*U;
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U) + Sm //+ Fe
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);

        K = 0.5*magSqr(U);
    }
