    MRF.correctBoundaryVelocity(U);

    // mixture.nu() = (1 + age*dimensionSet(0,0,-1,0,0))*mixture.nu();
    // mixture.nu() = 3*mixture.nu();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
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
    }
    

    // // FROM THE INTERNET
    // // for residential time
    // // if (ageExists)
    // {
    //     fvScalarMatrix ageEqn
    //     (
    //             fvm::ddt(age)
    //             + fvm::div(phi, age)
    //             - fvm::laplacian(turbulence->nuEff(), age)
    //             ==
    //             runTime.deltaT().value()
    //             // dimensionedScalar("ageSource", *age, runTime.deltaT().value())
    //             // dimensionedScalar("ageSource", *age.dimensions()*dimensionSet(0,0,-1,0,0), runTime.deltaT().value())
    //             //dimensionedScalar("ageSource", (*age).dimensions()*dimensionSet(0,0,-1,0,0), 1) // steady state
    //     );

    //     ageEqn.solve();
    // }



    // // WORKS: AGE, removed as it is calculated in the viscosity model
    // fvScalarMatrix ageEqn
    // (
    //         fvm::ddt(age)
    //         + fvm::div(phi, age)
    //         // - fvm::laplacian(turbulence->nuEff(), age)
    //         // + 5e-4
    //         // ==
    //         // runTime.deltaT().value()
    //         // dimensionedScalar("ageSource", age, runTime.deltaT().value())
    //         - dimensionedScalar(1)
    //         // + dimensionedScalar("ageSource", age.dimensions()*dimensionSet(0,0,-1,0,0), runTime.deltaT().value())
    //         //dimensionedScalar("ageSource", (*age).dimensions()*dimensionSet(0,0,-1,0,0), 1) // steady state
    // );

    // ageEqn.solve();
