/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UVTimeThickening.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(UVTimeThickening, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        UVTimeThickening,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::UVTimeThickening::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);


    const surfaceScalarField phiU(phi_); // /fvc::interpolate(deltaRho + deltaRho0));

    tmp<volScalarField> sr(strainRate());

    // dimensionedScalar dimCorr ("one",dimless/dimLength,1);
    volScalarField distFromCenter = sqrt(sqr(age_.mesh().C().component(1)) + 
                                         sqr(age_.mesh().C().component(2)));
    volScalarField xCutLower = (1 + sign(age_.mesh().C().component(0) - x1_)) / 2; //("one",dimless/sqr(dimLength),1)
    volScalarField xCutUpper = (1 + sign(x2_ - age_.mesh().C().component(0))) / 2; //("one",dimless/sqr(dimLength),1)
    
    // UV intensity field from the CT
    volScalarField intensityField = pow(d_, n_)/pow(distFromCenter, n_)*xCutLower*xCutUpper;

    // Info << phi_.dimensions() << endl;
    // Info << age_.dimensions() << endl;
    // Info << phiU.dimensions() << endl;    

    fvScalarMatrix ageEqn
    (
          fvm::ddt(age_)
        + fvm::div(phi_, age_)
        - intensityField
    );

    ageEqn.relax();
    ageEqn.solve();

    return
    (
        // min
        // (
        //     max
        //     (
        //         nu0_, 
        //         pow((1 + a_*age_),n_)*nu0_
        //     ),
        //     nuMax_
        // )
        min
        (
            max
            (
                nu0_, 
                nu0_ + (nuMax_ - nu0_)*(1 - exp(- a_ * age_))
            ),
            nuMax_
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::UVTimeThickening::UVTimeThickening
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    UVTimeThickeningCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    a_("a", dimless/dimTime, UVTimeThickeningCoeffs_),
    n_("n", dimless, UVTimeThickeningCoeffs_),
    d_("d", dimLength, UVTimeThickeningCoeffs_),
    x1_("x1", dimLength, UVTimeThickeningCoeffs_),
    x2_("x2", dimLength, UVTimeThickeningCoeffs_),
    nuMax_("nuMax", dimViscosity, UVTimeThickeningCoeffs_),
    nu0_("nu0", dimViscosity, UVTimeThickeningCoeffs_),
    
    // intensityField
    // (
    //     IOobject
    //     (
    //         "intensityField",
    //         U_.time().timeName(),
    //         U_.db(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     U_.mesh()
    // ),

    // const dictionary& myDict = mesh.controlDict();
    // scalar iWrite = myDict.lookupOrDefault<scalar>("printNu",0.);

    // If you use that, you cannot test if (iWrite == 1) as iWrite is a scalar so it might 
    // not be exactly equal to one (machine errors), just use a margin like (iWrite > 0.9 & iWrite < 1.1). 
    // Alternatively, you can try to replace scalar by int in the second line or by bool and use true/false instead of 0/1 but I am not sure if it works.

    age_
    (
        IOobject
        (
            "age",
            U_.time().timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    ),

    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::UVTimeThickening::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    UVTimeThickeningCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    UVTimeThickeningCoeffs_.lookup("a") >> a_;
    UVTimeThickeningCoeffs_.lookup("n") >> n_;
    UVTimeThickeningCoeffs_.lookup("d") >> d_;
    UVTimeThickeningCoeffs_.lookup("x1") >> x1_;
    UVTimeThickeningCoeffs_.lookup("x2") >> x2_;
    UVTimeThickeningCoeffs_.lookup("nuMax") >> nuMax_;
    UVTimeThickeningCoeffs_.lookup("nu0") >> nu0_;

    return true;
}


// ************************************************************************* //
