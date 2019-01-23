Info << "Reading field DVM Properties\n"
     << endl;
IOdictionary DVMProperties(
    IOobject(
        "DVMProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

Info << "Reading field rho\n"
     << endl;
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading field lambda\n"
     << endl;
volScalarField lambda(
    IOobject(
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volScalarField rhoOld(
    IOobject(
        "rhoOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("1", rho.dimensions(), 1));

volVectorField Uold(
    IOobject(
        "Uold",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedVector("1", U.dimensions(), vector(0, 0, 0)));

volScalarField lambdaOld(
    IOobject(
        "lambdaOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("1", lambda.dimensions(), 1));

volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("1", dimless / lambda.dimensions(), 1));

rhoOld = rho;
Uold = U;
lambdaOld = lambda;
T = 1.0 / lambda;

const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar maxCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

scalar convergeTol =
    runTime.controlDict().lookupOrDefault<scalar>("convergeTol", 1.0e-9);

label convergeCheckSteps =
    runTime.controlDict().lookupOrDefault("convergeCheckSteps", 1);

scalar lambdaChange = 1.0;
scalar rhoChange = 1.0;
scalar Uchange = 1.0;

fvDVM dvm(rho, U, lambda, DVMProperties);
