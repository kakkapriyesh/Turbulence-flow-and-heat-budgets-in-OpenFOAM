/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2010-08-02 Eelco van Vliet: 1st public version of wallHeatFluxIncompressible:
  http://www.cfd-online.com/Forums/openfoam-solving/66705-wallheatflux-bc-not-constant-after-restart.html#post269812

2012-05-21 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post362191
  «modified the standard wallHeatflux utility which comes default with OF into
  a version for incompressible flows. Also removed a bug out of the code.»

2012-06-26 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post368330
  «p is now not required anymore.»

2014-06-22: Bruno Santos: Adapted to OpenFOAM 2.2.x.

2018-06-15: Bruno Santos @ FSD blueCAPE Lda: Adapted to OpenFOAM 5.x.

-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    wallHeatFluxIncompressible

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.
    Based on wallHeatFlux with changes to allow it on incompressible flows
    Also removed a bug at the typeid checkline
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

#       include "createFields.H"
#       include "readTransportProperties.H"

         // update the turbulence fields
       turbulence->read();
}
      
{
       
    IOobject UPrime2MeanHeader
    (
        "UPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );


    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
          
     
  
    
    {
        Info<< "    Reading average field UPrime2Mean" << endl;
        const volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);

       // Info<< "    Reading average field kMean" << endl;
      //  const volScalarField kMean(kMeanHeader, mesh);

        Info<< "    Reading average field UMean" << endl;
        const volVectorField UMean(UMeanHeader, mesh);                          
     
        
        volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);                 


///////////////////////////////////////////////////////////////////////        
        
        
volScalarField TKEMean
        (
            IOobject
            (
                "TKEMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("TKEMean",UPrime2Mean.dimensions(),0)
        );
       

      
        
      
        
      


  volScalarField Ctt
        (
            IOobject
            (
                "Ctt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ctt",dimTemperature*dimTemperature/dimTime,0)
        );
     

        volScalarField Cut
        (
            IOobject
            (
                "Cut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Cut",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
         volScalarField Cvt
        (
            IOobject
            (
                "Cvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Cvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
         volScalarField Cwt
        (
            IOobject
            (
                "Cwt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Cwt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
  
 
volVectorField Cuug
(
    IOobject
    (
        "Cuug",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Cuug", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Cuu
(
    IOobject
    (
        "Cuu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Cuu", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);

volVectorField Cvvg
(
    IOobject
    (
        "Cvvg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Cvvg", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Cvv
(
    IOobject
    (
        "Cvv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Cvv", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);


volVectorField Cwwg
(
    IOobject
    (
        "Cwwg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Cwwg", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Cww
(
    IOobject
    (
        "Cww",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Cww", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);

volVectorField Cuvg
(
    IOobject
    (
        "Cuvg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Cuvg", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Cuv
(
    IOobject
    (
        "Cuv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Cuv", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);

volVectorField Cuwg
(
    IOobject
    (
        "Cuwg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Cuwg", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Cuw
(
    IOobject
    (
        "Cuw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Cuw", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);

volVectorField Cvwg
(
    IOobject
    (
        "Cvwg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Cvwg", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Cvw
(
    IOobject
    (
        "Cvw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Cvw", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);

volVectorField Ctkeg
(
    IOobject
    (
        "Ctkeg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedVector
      (
           "Ctkeg", 
           dimVelocity*dimVelocity/(dimLength),
vector (0,0,0)
      )
);

volScalarField Ctke
(
    IOobject
    (
        "Ctke",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
          
      dimensionedScalar
      (
           "Ctke", 
           dimVelocity*dimVelocity*dimVelocity/(dimLength),0
      )
);
  #include "createPhi.H"  
         
          singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );
        
        const volScalarField nuLam(turbulence->nu());                             

        Info<< "    Calculating mean turbulent kinetic energy TKEMean" << endl;
        TKEMean = 0.5 * (
                            UPrime2Mean.component(0)
                          + UPrime2Mean.component(3)
                          + UPrime2Mean.component(5)
                        );

     

    
        
    
        
             
        
     
  
  
     
    
  


//------------------------------------------ convection -----------------------------------------------------------------------

 Info<< "calculating convection" << endl;


Cuug = fvc::grad(UPrime2Mean.component(symmTensor::XX));
Cuu = Cuug.component(vector::X) * UMean.component(vector::X)+Cuug.component(vector::Y) * UMean.component(vector::Y)+Cuug.component(vector::Z) * UMean.component(vector::Z);

Cvvg = fvc::grad(UPrime2Mean.component(symmTensor::YY));
Cvv = Cvvg.component(vector::X) * UMean.component(vector::X)+Cvvg.component(vector::Y) * UMean.component(vector::Y)+Cvvg.component(vector::Z) * UMean.component(vector::Z);

Cwwg = fvc::grad(UPrime2Mean.component(symmTensor::ZZ));
Cww = Cwwg.component(vector::X) * UMean.component(vector::X)+Cwwg.component(vector::Y) * UMean.component(vector::Y)+Cwwg.component(vector::Z) * UMean.component(vector::Z);

Cuvg = fvc::grad(UPrime2Mean.component(symmTensor::XY));
Cuv = Cuvg.component(vector::X) * UMean.component(vector::X)+Cuvg.component(vector::Y) * UMean.component(vector::Y)+Cuvg.component(vector::Z) * UMean.component(vector::Z);

Cuwg = fvc::grad(UPrime2Mean.component(symmTensor::XZ));
Cuw = Cuwg.component(vector::X) * UMean.component(vector::X)+Cuwg.component(vector::Y) * UMean.component(vector::Y)+Cuwg.component(vector::Z) * UMean.component(vector::Z);

Cvwg = fvc::grad(UPrime2Mean.component(symmTensor::YZ));
Cvw = Cvwg.component(vector::X) * UMean.component(vector::X)+Cvwg.component(vector::Y) * UMean.component(vector::Y)+Cvwg.component(vector::Z) * UMean.component(vector::Z);

Ctkeg = fvc::grad(TKEMean);
Ctke = Ctkeg.component(vector::X) * UMean.component(vector::X)+Ctkeg.component(vector::Y) * UMean.component(vector::Y)+Ctkeg.component(vector::Z) * UMean.component(vector::Z);

  
     Cuu.write();
     Cvv.write();
     Cww.write();
     Cuv.write();
     Cuw.write();
     Cvw.write();
       Ctke.write();
     
   
    }
   
    
    Info<< "End" << endl;
}


// ************************************************************************* //
  
    return 0;
}

// ************************************************************************* //
