// 从twoPhaseProperties中读取各相导热系数等
IOdictionary Personset
(
    IOobject
    (
        "Personset",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
dimensionedScalar A_mush("A_mush", dimensionSet(1,-3,-1,0,0,0,0), Personset);
dimensionedScalar Tset("Tset", dimensionSet(0,0,0,1,0,0,0), Personset);
dimensionedScalar T0("T0", dimensionSet(0,0,0,1,0,0,0), 273.15);
dimensionedScalar T1("T1", dimensionSet(0,0,0,-1,0,0,0), 1);
//dimensionedScalar k_l("k_l", dimensionSet(0,0,0,-1,0,0,0), 
dimensionedScalar f0("f0", dimless, 0.044);
dimensionedScalar little("little", dimless, 0.001);
dimensionedScalar Pi("Pi", dimless, 3.1415926);
dimensionedScalar E("E", dimensionSet(1,1,-3,0,0,-1,0), 1);
dimensionedScalar Kb("Kb", dimless, 1.38e-7);
dimensionedScalar L("L", dimensionSet(0,2,-2,0,0,0,0), Personset);
dimensionedScalar nE("nE", dimensionSet(1,1,-3,0,0,-1,0), 1);
dimensionedScalar nS("nS", dimensionSet(0,0,-1,0,0,0,0), 1);
//dimensionedScalar IceCv("IceCv", dimensionSet(0,2,-2,-1,0,0,0), 2097);
dimensionedScalar IceCv("IceCv", dimensionSet(0,2,-2,-1,0,0,0), 2097);
dimensionedScalar IceKappa("IceKappa", dimensionSet(1,1,-3,-1,0,0,0), 2.1);
dimensionedScalar Icerho("Icerho", dimensionSet(1,-3,0,0,0,0,0), 917);
volScalarField& T = mixture.T();
volScalarField k_l
(
    IOobject
    (
        "k_l",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    8.0*Pi*0.0478*0.0478*0.0478*T0*T0*(2.0*T-Tset)*f0*exp((-4.7e-7)*mag(E)/nE)/(3.0*2.87*2.87*Kb*Tset*Tset*Tset*T*T*T1)
);
volScalarField f_ls
(
    IOobject
    (
        "f_ls",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.5 + 0.5*tanh(8.0*k_l*(T-T0))
);
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*(1-f_ls)
);
volScalarField SH
(
     IOobject
     (
         "SH",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
    (-4.0)*rho*k_l*Solid*L*(1.0-(tanh(8.0*k_l*(T-T0)))*(tanh(8.0*k_l*(T-T0))))*fvc::ddt(T)
);
volVectorField Sm
(
    IOobject
    (
        "Sm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*A_mush*(1.0-f_ls)/(f_ls*f_ls*f_ls+little)*U
);




