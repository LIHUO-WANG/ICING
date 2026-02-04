singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

Info<< "Reading nusselt parameters\n" << endl;
//wallTemp
dimensionedScalar Thot("Thot", dimTemperature, laminarTransport);
//charactic length
dimensionedScalar length("length", dimLength, laminarTransport);
//bottomWall
word heatWallName("down");
laminarTransport.readIfPresent("heatWallName", heatWallName);
