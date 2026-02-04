Info<< "Reading electric settings"<< nl  << endl;

IOdictionary electricProperties
(
    IOobject
    (
        "electricProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//int potentialNonorthogonalCorrectTime = electricProperties.subDict("settings").lookupOrDefault<int>("potentialNonorthogonalCorrectTime", 1);
//Info << "potentialNonorthogonalCorrectTime:" << potentialNonorthogonalCorrectTime << nl << endl;

//int electricFieldCorrectTime = electricProperties.subDict("settings").lookupOrDefault<int>("electricFieldCorrectTime", 1);
//Info << "electricFieldCorrectTime:" << electricFieldCorrectTime << nl << endl;

//electric charge
dimensionedScalar e0("e0", dimensionSet(0,0,1,0,0,1,0), 1.6e-19);
//Boltzmann  constant
dimensionedScalar B("B", dimensionSet(1,2,-2,-1,0,0,0), 1.38e-23);
//ion mobility
dimensionedScalar Kp("Kp", dimensionSet(-1,0,2,0,0,1,0), electricProperties);
dimensionedScalar Kn("Kn", dimensionSet(-1,0,2,0,0,1,0), electricProperties);
dimensionedScalar Dp("Dp", dimensionSet(0,2,-1,0,0,0,0), electricProperties);
dimensionedScalar Dn("Dn", dimensionSet(0,2,-1,0,0,0,0), electricProperties);
//pemitivity
dimensionedScalar eps_Air("eps_Air", dimensionSet(-1,-3,4,0,0,2,0), electricProperties);
dimensionedScalar eps_Water("eps_Water", dimensionSet(-1,-3,4,0,0,2,0), electricProperties);
dimensionedScalar eps0("eps0", dimensionSet(-1,-3,4,0,0,2,0), 8.85e-12);
//conductivity
dimensionedScalar sig("sig", dimensionSet(-1,-3,3,0,0,2,0), electricProperties);
//density
//dimensionedScalar rho("rho", dimensionSet(1,-3,0,0,0,0,0), electricProperties);
//ion concentration without electric field
dimensionedScalar rho0(sig/(Kp+Kn));

//recombination coffecient
//dimensionedScalar eps("eps", dimensionSet(-1,-3,4,0,0,2,0));
//eps = 1/(alpha1/eps_Water + (1-alpha1)/eps_Air);

//dimensionedScalar krx((Kp+Kn)/eps);
