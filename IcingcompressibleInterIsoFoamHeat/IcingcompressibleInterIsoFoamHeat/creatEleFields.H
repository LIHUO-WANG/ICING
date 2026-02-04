#include"readEleDict.H"

Info<< "Creating electric Fields\n" << endl;

//potential
Info<< "Reading field Ep\n" << endl;
volScalarField Ep
(
    IOobject
    (
        "Ep",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//npchageDensity  
Info<< "Reading field rhop\n" << endl;
volScalarField rhop
(
   IOobject
   (
       "rhop",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
   ),
   mesh
);

//nnchageDensity  
Info<< "Reading field rhon\n" << endl;
volScalarField rhon
(
    IOobject
    (
        "rhon",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rhoe
(
   IOobject
   (
       "rhoe",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
   ),
   rhop-rhon
);

//electric field flux  
Info<< "Calculating field rhoeFlux\n" << endl;
surfaceScalarField rhoeFluxP
(
   IOobject
   (
       "rhoeFluxP",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::NO_WRITE
   ),
   Kp*mesh.magSf()*fvc::snGrad(Ep)
);

surfaceScalarField rhoeFluxN
(
    IOobject
    (
        "rhoeFluxN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Kn*mesh.magSf()*fvc::snGrad(Ep)
);

//electric field 
Info<< "Calculating field E\n" << endl;
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(Ep)
);
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1/(alpha1/eps_Water + (1-alpha1)/eps_Air)
);
dimensionedScalar Euni("Euni", dimensionSet(1,1,-3,0,0,-1,0), 1);
volScalarField krx
(
    IOobject
    (
        "krx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Kp+Kn)/eps
);
volVectorField Dis
(
    IOobject
    (
        "Dis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   E*eps,
   E.boundaryField().types()//因为我们后期是重建的Dis场，所以Dis的边界条件应该是通过内部的场外推边界的值，通过浏览代码可以发现，grad命令得到的场其边界条件就是外推，因此这里直接进行复制。
);

//electric force
Info<< "Caluclating field Fe\n" << endl;

volVectorField Fe
(
    IOobject
    (
        "Fe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    // rhon*E
    rhoe*E
);

volScalarField uniField
(
    IOobject
    (
        "uniField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "interfaceCellMaker",
         dimless,
         1
    )
);

volScalarField Efunction
(
    IOobject
    (
        "Efunction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "interfaceCellMaker",
         dimless,
         1
    )
);

volScalarField sourceTerm
(
   IOobject
   (
       "sourceTerm",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::NO_WRITE
   ),
   Efunction*krx*rho0*rho0-krx*rhop*rhon
);


