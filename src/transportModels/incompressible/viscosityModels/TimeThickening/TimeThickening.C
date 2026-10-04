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

#include "TimeThickening.H"
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
    defineTypeNameAndDebug(TimeThickening, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        TimeThickening,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::TimeThickening::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);


    const surfaceScalarField phiU(phi_); // /fvc::interpolate(deltaRho + deltaRho0));

    tmp<volScalarField> sr(strainRate());

    // Info << phi_.dimensions() << endl;
    // Info << age_.dimensions() << endl;
    // Info << phiU.dimensions() << endl;    

    fvScalarMatrix ageEqn
    (
          fvm::ddt(age_)
        + fvm::div(phi_, age_)
        - dimensionedScalar("one",dimless,1)
    );

    ageEqn.relax();
    ageEqn.solve();

    return
    (
        min
        (
            max
            (
                nu0_, 
                pow((1 + a_*age_),n_)*nu0_
            ),
            nuMax_
        )
        
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::TimeThickening::TimeThickening
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    TimeThickeningCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    a_("a", dimless/dimTime, TimeThickeningCoeffs_),
    n_("n", dimless, TimeThickeningCoeffs_),
    nuMax_("nuMax", dimViscosity, TimeThickeningCoeffs_),
    nu0_("nu0", dimViscosity, TimeThickeningCoeffs_),
    age_
    (
        IOobject
        (
            "age",
            U_.time().timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::TimeThickening::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    TimeThickeningCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    TimeThickeningCoeffs_.lookup("a") >> a_;
    TimeThickeningCoeffs_.lookup("n") >> n_;
    TimeThickeningCoeffs_.lookup("nuMax") >> nuMax_;
    TimeThickeningCoeffs_.lookup("nu0") >> nu0_;

    return true;
}


// ************************************************************************* //
