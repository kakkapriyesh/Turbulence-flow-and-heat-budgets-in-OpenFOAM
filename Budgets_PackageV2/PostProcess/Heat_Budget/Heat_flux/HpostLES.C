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

    // IOobject kMeanHeader
   // (
   //     "kMean",
   //     runTime.timeName(),
   //     mesh,
   //     IOobject::MUST_READ
   // );

    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
          
   
   // if
   // (
  //      UPrime2MeanHeader.headerOk()
  //   && kMeanHeader.headerOk()
 //    && UMeanHeader.headerOk()
  //   && turbDiffMeanHeader.headerOk()
 //    && SGSDiffMeanHeader.headerOk()     
 //   )
    
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
            dimensionedScalar("dimkMean",UPrime2Mean.dimensions(),0)
        );

        volScalarField resLES
        (
            IOobject
            (
                "resLES",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimless",dimless,0)
        );

        volScalarField TKEMeanProd
        (
            IOobject
            (
                "TKEMeanProd",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimTKEMeanProd",UPrime2Mean.dimensions()/dimTime,0)
        );
        
       
        volVectorField Gt
(
    IOobject
    (
        "Gt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Gt",dimTemperature/dimLength,vector (0,0,0))
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
     
         volScalarField Ptt
        (
            IOobject
            (
                "Ptt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ptt",dimTemperature*dimTemperature/dimTime,0)
        );
          volScalarField Dtt
        (
            IOobject
            (
                "Dtt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dtt",dimTemperature*dimTemperature/(dimLength*dimLength),0)
        );
           volScalarField Ttt
        (
            IOobject
            (
                "Ttt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ttt",dimTemperature*dimTemperature/(dimTime),0)
        );
           volScalarField Mtt
        (
            IOobject
            (
                "Mtt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Mtt",dimTemperature*dimTemperature/(dimLength*dimLength),0)
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
            dimensionedScalar("Cut",dimTemperature*(dimLength/(dimTime*dimTime)),0)
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
        volScalarField Put
        (
            IOobject
            (
                "Put",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Put",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
 volScalarField Pvt
        (
            IOobject
            (
                "Pvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Pvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
 volScalarField Pwt
        (
            IOobject
            (
                "Pwt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Pwt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
          volScalarField Dvt
        (
            IOobject
            (
                "Dvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
           volScalarField Dut
        (
            IOobject
            (
                "Dut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dut",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
volScalarField Dwt
(
            IOobject
            (
                "Dwt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dwt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
          volScalarField Tvt
        (
            IOobject
            (
                "Tvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
        volScalarField viscDiffusionMean
        (
            IOobject
            (
                "viscDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("viscDiffusionMean",UPrime2Mean.dimensions()/dimTime,0)
        ); 
          Info<< "Reading field Tmean\n" << endl;
    volScalarField TMean
    (
        IOobject
        (
            "TMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
     Info<< "Reading field TUmean\n" << endl;
    volVectorField TUMean
    (
        IOobject
        (
            "TUMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
      Info<< "Reading field TTUmean\n" << endl;
    volVectorField TTUMean
    (
        IOobject
        (
            "TTUMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
  Info<< "Reading field TTDissMean\n" << endl;
    volVectorField TTDissMean
    (
        IOobject
        (
            "TTDissMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
  
  
     Info<< "Reading field UUTMean\n" << endl;
    volScalarField WWUMean
    (
        IOobject
        (
            "WWUMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
     volScalarField WWVMean
    (
        IOobject
        (
            "WWVMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

volVectorField GradTUPMean
(
    IOobject
    (
        "GradTUPMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
           
);

volTensorField UUTMean
(
    IOobject
    (
        "UUTMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
           
);

 volVectorField uuT
        (
            IOobject
            (
                "uuT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("uuT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );
volVectorField uvT
        (
            IOobject
            (
                "uvT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("uvT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );
volVectorField uwT
        (
            IOobject
            (
                "uwT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("uwT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );
volVectorField vuT
        (
            IOobject
            (
                "vuT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("vuT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );
volVectorField vvT
        (
            IOobject
            (
                "vvT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("vvT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );
volVectorField vwT
        (
            IOobject
            (
                "vwT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("vwT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );

volVectorField wuT
        (
            IOobject
            (
                "wuT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("wuT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );

volVectorField wvT
        (
            IOobject
            (
                "wvT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("wvT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );

        

volVectorField wwT
        (
            IOobject
            (
                "wwT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
        
    dimensionedVector("wwT",(dimTemperature*dimLength)/(dimTime*dimTime),vector (0,0,0))
        );

volScalarField TDuT
(
    IOobject
    (
        "TDuT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
 dimensionedScalar("TDuT",(dimTemperature*dimLength)/(dimTime*dimTime),0)
           
);
volScalarField TDvT
(
    IOobject
    (
        "TDvT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
 dimensionedScalar("TDvT",(dimTemperature*dimLength)/(dimTime*dimTime),0)
           
);
volScalarField TDwT
(
    IOobject
    (
        "TDwT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
 dimensionedScalar("TDwT",(dimTemperature*dimLength)/(dimTime*dimTime),0)
           
);

volScalarField TPDu
(
    IOobject
    (
        "TPDu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
 dimensionedScalar("TPDu",(dimTemperature*dimLength)/(dimTime*dimTime),0)
           
);

volScalarField TPDv
(
    IOobject
    (
        "TPDv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
 dimensionedScalar("TPDv",(dimTemperature*dimLength)/(dimTime*dimTime),0)
           
);

volScalarField TPDw
(
    IOobject
    (
        "TPDw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedScalar("TPDv",(dimTemperature*dimLength)/(dimTime*dimTime),0)
  
           
);

volVectorField TGradp
( 
    IOobject
    (
        "TGradp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
),
    mesh

 
);


     volScalarField WWWMean
    (
        IOobject
        (
            "WWWMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

        


volTensorField MD1Mean
(
    IOobject
    (
        "MD1Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
           
);

volTensorField MD2Mean
(
    IOobject
    (
        "MD2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
           
);



volVectorField MD1Tux
(
    IOobject
    (
        "MD1Tux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD1Tux",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);

volVectorField MD1Tuy
(
    IOobject
    (
        "MD1Tuy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD1Tuy",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);
volVectorField MD1Tuz
(
    IOobject
    (
        "MD1Tuz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
   dimensionedVector("MD1Tuz",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);

volVectorField MD1Tvx
(
    IOobject
    (
        "MD1Tvx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD1Tvx",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);

volVectorField MD1Tvy
(
    IOobject
    (
        "MD1Tvy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD1Tvy",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);
volVectorField MD1Tvz
(
    IOobject
    (
        "MD1Tvz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
   dimensionedVector("MD1Tvz",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);
volVectorField MD1Twx
(
    IOobject
    (
        "MD1Twx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD1Twx",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);

volVectorField MD1Twy
(
    IOobject
    (
        "MD1Twy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD1Twy",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);
volVectorField MD1Twz
(
    IOobject
    (
        "MD1Twz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
   dimensionedVector("MD1Twz",((dimTemperature)/(dimTime*dimLength)),vector (0,0,0))
      
);


volVectorField MD2utx
(
    IOobject
    (
        "MD2utx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2utx",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
      
);

volVectorField MD2uty
(
    IOobject
    (
        "MD2uty",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2uty",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
      
);

volVectorField MD2utz
(
    IOobject
    (
        "MD2utz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2utz",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
       
);

volVectorField MD2wtx
(
    IOobject
    (
        "MD2wtx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2wtx",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
      
);

volVectorField MD2wty
(
    IOobject
    (
        "MD2wty",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2wty",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
      
);

volVectorField MD2wtz
(
    IOobject
    (
        "MD2wtz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2wtz",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
       
);

volVectorField MD2vtx
(
    IOobject
    (
        "MD2vtx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2vtx",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
      
);

volVectorField MD2vty
(
    IOobject
    (
        "MD2vty",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2vty",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
      
);

volVectorField MD2vtz
(
    IOobject
    (
        "MD2vtz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
    dimensionedVector("MD2vtz",(dimTemperature)/(dimTime*dimLength),vector (0,0,0))
       
);
volScalarField MDu1
(
    IOobject
    (
        "MDu1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
dimensionedScalar("MDu1",(dimTemperature*dimLength)/(dimTime*dimTime),0)
      
);
volScalarField MDv1
(
    IOobject
    (
        "MDv1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
dimensionedScalar("MDv1",(dimTemperature*dimLength)/(dimTime*dimTime),0)
      
);
volScalarField MDw1
(
    IOobject
    (
        "MDw1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
dimensionedScalar("MDw1",(dimTemperature*dimLength)/(dimTime*dimTime),0)
      
);
volScalarField MDu2
(
    IOobject
    (
        "MDu2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
dimensionedScalar("MDu2",(dimTemperature*dimLength)/(dimTime*dimTime),0)
      
);
volScalarField MDv2
(
    IOobject
    (
        "MDv2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
dimensionedScalar("MDv2",(dimTemperature*dimLength)/(dimTime*dimTime),0)
      
);
volScalarField MDw2
(
    IOobject
    (
        "MDw2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
        
dimensionedScalar("MDw2",(dimTemperature*dimLength)/(dimTime*dimTime),0)
      
);
volScalarField nutMean
    (
        IOobject
        (
            "nutMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

volScalarField alphatMean
    (
        IOobject
        (
            "alphatMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

  #include "createPhi.H"  
        //~ #include "createFields.H" 
        
          singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );
        
        const volScalarField nuLam(turbulence->nu());             
volScalarField alpha("alpha", turbulence->nu()/0.71);  //0.71 is Prandtl number for air               

        Info<< "    Calculating mean turbulent kinetic energy TKEMean" << endl;
        TKEMean = 0.5 * (
                            UPrime2Mean.component(0)
                          + UPrime2Mean.component(3)
                          + UPrime2Mean.component(5)
                        );

     

        Info<< "    Calculating Mean TKE production TKEMeanProd" << endl;
    
        TKEMeanProd = - UPrime2Mean && fvc::grad(UMean);
        
        
        
        Info<< "    Calculating viscDiffusionMean" << endl;
        viscDiffusionMean = nuLam*fvc::laplacian(TKEMean);         
        
     
        
        
        Info<< "    Calculating convection" << endl;
        Gt= fvc::grad(TMean);
    
        Info<< "    Calculating convection ut" << endl;
        Cut= UMean & fvc::grad(TUMean.component(0));
        Cvt= UMean & fvc::grad(TUMean.component(1));
        Cwt= UMean & fvc::grad(TUMean.component(2));
        Info<< "    Calculating Production ut" << endl;
        Put= -((TUMean & fvc::grad(UMean.component(0)))+(UPrime2Mean.component(0)*Gt.component(0)+UPrime2Mean.component(1)*Gt.component(1)+UPrime2Mean.component(2)*Gt.component(2)));
        Pvt= -((TUMean & fvc::grad(UMean.component(1)))+(UPrime2Mean.component(1)*Gt.component(0)+UPrime2Mean.component(3)*Gt.component(1)+UPrime2Mean.component(4)*Gt.component(2)));
        Pwt= -((TUMean & fvc::grad(UMean.component(2)))+(UPrime2Mean.component(2)*Gt.component(0)+UPrime2Mean.component(4)*Gt.component(1)+UPrime2Mean.component(5)*Gt.component(2)));
       Info<< "    Calculating Dissipation ut" << endl;
   Dut= -(nutMean+nuLam+alpha+alphatMean)*GradTUPMean.component(0);
   Dvt= -(nutMean+nuLam+alpha+alphatMean)*GradTUPMean.component(1);
   Dwt= -(nutMean+nuLam+alpha+alphatMean)*GradTUPMean.component(2);
        Info<< "    Calculating T diffusion ut" << endl;
   uuT= fvc::grad(UUTMean.component(0)) ; 
   uvT= fvc::grad(UUTMean.component(1)) ;
   uwT= fvc::grad(UUTMean.component(2)) ; 
   vuT= fvc::grad(UUTMean.component(3)) ; 
   vvT= fvc::grad(UUTMean.component(4)) ;
   vwT= fvc::grad(UUTMean.component(5)) ; 
   wuT= fvc::grad(UUTMean.component(6)) ; 
   wvT= fvc::grad(UUTMean.component(7)) ;
   wwT= fvc::grad(UUTMean.component(8)) ; 

   TDuT= -(uuT.component(0)+uvT.component(1)+uwT.component(2));
   TDvT= -(vuT.component(0)+vvT.component(1)+vwT.component(2));
   TDwT= -(wuT.component(0)+wvT.component(1)+wwT.component(2));


        Info<< "    Calculating Temperature pressure diffusion ut" << endl;
TPDu = TGradp.component(0);
TPDv =  TGradp.component(1);
TPDw =  TGradp.component(2);

 Info<< "    Calculating Molecular diffusion" << endl;

 MD1Tux= fvc::grad(MD1Mean.component(0)) ; 
 MD1Tvx= fvc::grad(MD1Mean.component(1)) ;
 MD1Twx= fvc::grad(MD1Mean.component(2)) ;
 MD1Tuy= fvc::grad(MD1Mean.component(3)) ; 
 MD1Tvy= fvc::grad(MD1Mean.component(4)) ;
 MD1Twy= fvc::grad(MD1Mean.component(5)) ;
 MD1Tuz= fvc::grad(MD1Mean.component(6)) ; 
 MD1Tvz= fvc::grad(MD1Mean.component(7)) ;
 MD1Twz= fvc::grad(MD1Mean.component(8)) ;

 MD2utx= fvc::grad(MD2Mean.component(0)) ; 
 MD2uty= fvc::grad(MD2Mean.component(1)) ;
 MD2utz= fvc::grad(MD2Mean.component(2)) ; 
 MD2vtx= fvc::grad(MD2Mean.component(3)) ;
 MD2vty= fvc::grad(MD2Mean.component(4)) ;
 MD2vtz= fvc::grad(MD2Mean.component(5)) ;
 MD2wtx= fvc::grad(MD2Mean.component(6)) ; 
 MD2wty= fvc::grad(MD2Mean.component(7)) ;
 MD2wtz= fvc::grad(MD2Mean.component(8)) ;

MDu1 = (nuLam+nutMean)*(MD1Tux.component(0)+MD1Tuy.component(1)+MD1Tuz.component(2));
MDu2=(alpha+alphatMean)*(MD2utx.component(0)+MD2uty.component(1)+MD2utz.component(2));
MDv1 = (nuLam+nutMean)*(MD1Tvx.component(0)+MD1Tvy.component(1)+MD1Tvz.component(2));
MDv2=(alpha+alphatMean)*(MD2vtx.component(0)+MD2vty.component(1)+MD2vtz.component(2));
MDw1 = (nuLam+nutMean)*(MD1Twx.component(0)+MD1Twy.component(1)+MD1Twz.component(2));
MDw2=(alpha+alphatMean)*(MD2wtx.component(0)+MD2wty.component(1)+MD2wtz.component(2)); // doubt in constant coefficient
  
       
      Cut.write();
      Cvt.write();
      Cwt.write();
      Put.write();
      Pvt.write();
      Pwt.write();
      Dut.write();
      Dvt.write();
      Dwt.write();
      TDuT.write();
      TDvT.write();
      TDwT.write();
      TPDu.write();
      TPDv.write();
      TPDw.write();
      MDu1.write();
	 MDv1.write();
	 MDw1.write();
      MDu2.write();
	 MDv2.write();
	 MDw2.write();


   
    }
   
    
    Info<< "End" << endl;
}


// ************************************************************************* //
  

    return 0;
}

// ************************************************************************* //
