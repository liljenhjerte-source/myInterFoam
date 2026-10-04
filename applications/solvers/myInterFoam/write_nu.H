{
        volScalarField shearRate
        (
                IOobject
                (
                        "shearRate",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                mag(fvc::grad(U))
        );
        shearRate.write();

        volScalarField strainRate
        (
                IOobject
                (
                        "strainRate",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                Foam::sqrt(2.0)*mag(symm(fvc::grad(U)))
        );
        strainRate.write();

    volScalarField nu
    (
            IOobject
            (
                    "nu",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
            ),
            // interfaceProperties.nu()
            // tmp<volScalarField> nu()
            mixture.nu()
    ); 
    nu.write();
}