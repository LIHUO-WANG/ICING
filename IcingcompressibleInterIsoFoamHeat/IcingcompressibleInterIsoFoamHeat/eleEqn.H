
 // solving the Poisson Equation
 fvScalarMatrix EpEqn
 (
     fvm::laplacian(eps,Ep) + rhoe
     //fvm::laplacian(eps,Ep) + rhon
 );
 EpEqn.solve();
 
//#include "reconstructEle.H"				//Calculation of Electric Field
Ep.correctBoundaryConditions();
E=-fvc::grad(Ep);
E.correctBoundaryConditions();

rhoeFluxP = phi-Kp*mesh.magSf()*fvc::snGrad(Ep);  	//Definition of charge flux
rhoeFluxN = phi+Kn*mesh.magSf()*fvc::snGrad(Ep);  	//Definition of charge flux

//rhoeFluxP = -alphaSuf*Kp*mesh.magSf()*fvc::snGrad(Ep);  	//Definition of charge flux
//rhoeFluxN = alphaSuf*Kn*mesh.magSf()*fvc::snGrad(Ep);  	//Definition of charge flux

//rhoeFluxP=phiCorrected;
//rhoeFluxN=phiCorrected;
krx = (Kp+Kn)/eps;
Efunction=uniField+(e0*e0*e0*mag(E))/(8.0*3.14*eps*B*B*T*T);

//sourceTerm=alpha1*Efunction*krx*rho0*rho0-krx*rhop*rhon;
sourceTerm=Efunction*krx*rho0*rho0-krx*rhop*rhon;

fvScalarMatrix rhopEqn
(
   fvm::ddt(rhop)   
   + fvm::div(rhoeFluxP, rhop)	
   - fvm::laplacian(Dp, rhop)
   //== sourceTerm-fvm::Sp(krx*rhon,rhop)
   ==sourceTerm
);
rhopEqn.solve();
rhop.correctBoundaryConditions();

// solving the positive charge conservation equation
fvScalarMatrix rhonEqn
(
    fvm::ddt(rhon)
    + fvm::div(rhoeFluxN, rhon)			
    - fvm::laplacian(Dn, rhon)			
    //== sourceTerm-fvm::Sp(krx*rhop,rhon)
    ==sourceTerm
);
rhonEqn.solve();
rhon.correctBoundaryConditions();

rhoe=rhop-rhon;
rhoe.correctBoundaryConditions();
//Fe=rhoe*E;
eps = 1/(alpha1/eps_Water + (1-alpha1)/eps_Air);
Fe=rhoe*E+0.5*((E&E)*fvc::grad(eps));
