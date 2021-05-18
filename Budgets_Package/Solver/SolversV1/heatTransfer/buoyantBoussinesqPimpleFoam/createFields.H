Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field TMean\n" << endl;
volScalarField TMean
(
    IOobject
    (
        "TMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field UMean\n" << endl;
volVectorField UMean
(
    IOobject
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field pMean\n" << endl;
volScalarField pMean
(
    IOobject
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField UPrime
(
    IOobject
    (
        "UPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedVector
      (
           "UPrime", 
           dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volVectorField PU
(
    IOobject
    (
        "PU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedVector
      (
           "PU", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
//~ 
volScalarField pPrime
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pPrime", 
           dimVelocity*dimVelocity, // Or Umean.dimensions()
           scalar(0.0)
      )    
);
volScalarField TPrime
(
    IOobject
    (
        "TPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "TPrime", 
           dimTemperature, // Or Umean.dimensions()
           scalar(0.0)
      )    
);
//~ 

volSymmTensorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "B", 
           dimVelocity*dimVelocity,
           symmTensor::zero
      )
);


volSymmTensorField SGSstrainTensor
(
    IOobject
    (
        "SGSstrainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "SGSstrainTensor", 
           dimVelocity/dimLength,
           symmTensor::zero
      )
);

volScalarField viscDiss
(
    IOobject
    (
        "viscDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "viscDiss", 
           dimVelocity*dimVelocity*dimVelocity/dimLength,
           0.0
      )
);

volScalarField pressDiff
(
    IOobject
    (
        "pressDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pressDiff", 
           dimVelocity*dimLength/(dimTime*dimTime),
           0.0
      )
);


volVectorField turbDiff
(
   IOobject
    (
        "turbDiff",
       runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
          "turbDiff", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volScalarField SGSDiss
(
    IOobject
    (
        "SGSDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "SGSDiss", 
           dimVelocity*dimVelocity*dimVelocity/dimLength,
           0.0
      )
);


volVectorField SGSDiff
(
    IOobject
    (
        "SGSDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "SGSDiff", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volSymmTensorField DP
(
    IOobject
    (
        "DP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "DP", 
           dimVelocity/(dimVelocity*dimTime), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField TijP
(
    IOobject
    (
        "TijP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "TijP", 
           (dimVelocity*dimVelocity), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdiffuTP
(
    IOobject
    (
        "SGSdiffuTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdiffuTP", 
           (dimVelocity*dimVelocity*dimVelocity), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdiffvTP
(
    IOobject
    (
        "SGSdiffvTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdiffvTP", 
           (dimVelocity*dimVelocity*dimVelocity), // Or Umean.dimensions()
           symmTensor::zero
      )
);

volSymmTensorField SGSdiffwTP
(
    IOobject
    (
        "SGSdiffwTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdiffwTP", 
           (dimVelocity*dimVelocity*dimVelocity), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdissuxTP
(
    IOobject
    (
        "SGSdissuxTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdissuxTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdissvxTP
(
    IOobject
    (
        "SGSdissvxTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdissvxTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdisswxTP
(
    IOobject
    (
        "SGSdisswxTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdisswxTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdissuyTP
(
    IOobject
    (
        "SGSdissuyTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdissuyTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdissvyTP
(
    IOobject
    (
        "SGSdissvyTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdissvyTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdisswyTP
(
    IOobject
    (
        "SGSdisswyTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdisswyTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdissuzTP
(
    IOobject
    (
        "SGSdissuzTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdissuzTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdissvzTP
(
    IOobject
    (
        "SGSdissvzTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdissvzTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
volSymmTensorField SGSdisswzTP
(
    IOobject
    (
        "SGSdisswzTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedSymmTensor
      (
           "SGSdisswzTP", 
           (dimVelocity*dimVelocity*dimVelocity/(dimLength)), // Or Umean.dimensions()
           symmTensor::zero
      )
);
/*volVectorField X
(
    IOobject
    (
        "X",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "X", 
           dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);*/
#include "createPhi.H"

#include "readTransportProperties.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - TRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volSymmTensorField PS
(
    IOobject
    (
        "PS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "PS", 
           dimVelocity*dimLength/(dimTime*dimTime),
           symmTensor::zero
      )
);

/*volTensorField turbDiffTest
(
    IOobject
    (
        "turbDiffTest",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedTensor
      (
           "turbDiffTest", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);*/


volSymmTensorField PD
(
    IOobject
    (
        "PD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "PD", 
           dimVelocity*dimLength/(dimTime*dimTime),
           symmTensor::zero
      )
);
volScalarField UUU
(
    IOobject
    (
        "UUU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "UUU", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField UUV
(
    IOobject
    (
        "UUV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "UUV", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField UUW
(
    IOobject
    (
        "UUW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "UUW", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField VVV
(
    IOobject
    (
        "VVV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "VVV", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField VVU
(
    IOobject
    (
        "VVU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "VVU", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField VVW
(
    IOobject
    (
        "VVW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "VVW", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField WWW
(
    IOobject
    (
        "WWW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "WWW", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField WWV
(
    IOobject
    (
        "WWV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "WWV", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField WWU
(
    IOobject
    (
        "WWU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "WWU", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);
volScalarField UVW
(
    IOobject
    (
        "UVW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "UVW", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           0.0
      )
);

volTensorField GradUP
(
    IOobject
    (
        "GradUP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "GradUP", 
           dimVelocity/dimLength,
           tensor::zero
      )
);
volTensorField tkeDiss
(
    IOobject
    (
        "tkeDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "tkeDiss", 
           (dimVelocity*dimVelocity/dimTime),
           //(dimVelocity/(dimLength*dimTime)),
           tensor::zero
           //0.0
      )
);


volScalarField dissUV
(
    IOobject
    (
        "dissUV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "dissUV", 
           (dimVelocity*dimVelocity/dimTime),
           0.0
      )
);
volScalarField dissUW
(
    IOobject
    (
        "dissUW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "dissUW", 
           (dimVelocity*dimVelocity/dimTime),
           0.0
      )
);
volScalarField dissVW
(
    IOobject
    (
        "dissVW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "dissVW", 
           (dimVelocity*dimVelocity/dimTime),
           0.0
      )
);
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);
// Temp variables.................................................

volVectorField TU
(
    IOobject
    (
        "TU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "TU", 
           dimVelocity*dimTemperature, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volVectorField TTU
(
    IOobject
    (
        "TTU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "TTU", 
           dimVelocity*dimTemperature*dimTemperature, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volVectorField GradTP
(
    IOobject
    (
        "GradTP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "GradTP", 
           dimTemperature/dimLength, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volVectorField GradTUP
(
    IOobject
    (
        "GradTUP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "GradTUP", 
           dimVelocity*dimTemperature/(dimLength*dimLength), // Or Umean.dimensions()
           vector (0,0,0)
      )
);


volVectorField GradTUPL
(
    IOobject
    (
        "GradTUPL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "GradTUPL", 
           dimVelocity*dimTemperature/(dimLength*dimLength), // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volVectorField GradTUPL2
(
    IOobject
    (
        "GradTUPL2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "GradTUPL2", 
           dimVelocity*dimTemperature/(dimLength*dimLength), // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volScalarField GradTU
(
    IOobject
    (
        "GradTU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "GradTU", 
         dimVelocity*dimTemperature/(dimLength*dimLength),0 // Or Umean.dimensions()
           
      )
);
volVectorField GradTUP1
(
    IOobject
    (
        "GradTUP1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "GradTUP1", 
           dimVelocity*dimTemperature/(dimLength*dimLength), // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volScalarField GradTV
(
    IOobject
    (
        "GradTV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "GradTV", 
         dimVelocity*dimTemperature/(dimLength*dimLength),0 // Or Umean.dimensions()
           
      )
);
volScalarField GradTW
(
    IOobject
    (
        "GradTW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "GradTW", 
         dimVelocity*dimTemperature/(dimLength*dimLength),0 // Or Umean.dimensions()
           
      )
);
volVectorField TTDiss
(
    IOobject
    (
        "TTDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "TTDiss", 
           (dimTemperature/dimLength)*(dimTemperature/dimLength), // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volVectorField TTDissL
(
    IOobject
    (
        "TTDissL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "TTDissL", 
           (dimTemperature/dimLength)*(dimTemperature/dimLength), // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volVectorField TTDiffL
(
    IOobject
    (
        "TTDiffL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "TTDiffL", 
           (dimTemperature*dimTemperature/(dimLength)), // Or Umean.dimensions()
           vector (0,0,0)
      )
);


volTensorField UUT
(
    IOobject
    (
        "UUT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedTensor
      (
           "UUT", 
           dimVelocity*dimVelocity*dimTemperature, // Or Umean.dimensions()
           tensor::zero
      )
);
volVectorField TGradp
(
    IOobject
    (
        "TGradp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "TGradp", 
           (dimTemperature*dimLength)/(dimTime*dimTime), // Or Umean.dimensions()
           vector (0,0,0)
      )
);
volVectorField Gradp
(
    IOobject
    (
        "Gradp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "Gradp", 
           (dimLength)/(dimTime*dimTime), // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volTensorField MD1
(
    IOobject
    (
        "MD1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "MD1", 
           dimTemperature/dimTime,
           tensor::zero
      )
);

volTensorField MD1L
(
    IOobject
    (
        "MD1L",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "MD1L", 
           dimTemperature/dimTime,
           tensor::zero
      )
);







volTensorField MD2
(
    IOobject
    (
        "MD2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "MD2", 
           dimTemperature/dimTime,
           tensor::zero
      )
);


volTensorField MD2L
(
    IOobject
    (
        "MD2L",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "MD2L", 
           dimTemperature/dimTime,
           tensor::zero
      )
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
