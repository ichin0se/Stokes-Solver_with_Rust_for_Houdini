#include "SIM_Stokes.hpp"
#include "./util/eigen.h"
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_PerfMonAutoEvent.h>
#include <UT/UT_VoxelArray.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_FieldSampler.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_RawIndexField.h>
#include <SIM/SIM_Object.h>
#include <GAS/GAS_SubSolver.h>
#include <CE/CE_Vector.h>
#include <CE/CE_SparseMatrix.h>

using namespace Stokes;

void initializeSIM(void *) {
    Eigen::setNbThreads(12);
    #ifndef NDEBUG
    if ( Eigen::nbThreads() > 1 )
        std::cout << "Eigen is multithreaded" << std::endl;
    #endif
    IMPLEMENT_DATAFACTORY(SIM_Stokes);
}

SIM_Stokes::SIM_Stokes(const SIM_DataFactory *factory) : BaseClass(factory) {
}

SIM_Stokes::~SIM_Stokes() {
}

const SIM_DopDescription * SIM_Stokes::getDopDescription() {
    static PRM_Name theVelocityName(GAS_NAME_VELOCITY, "Velocity Field");
    static PRM_Default theVelocityDefault(0, "vel");
    static PRM_Name theViscosityName("viscosity", "Viscosity Field");
    static PRM_Default theViscosityDefault(0, "viscosity");
    static PRM_Name theCollisionName(GAS_NAME_COLLISION, "Collision Field");
    static PRM_Default theCollisionDefault(0, "collision");
    static PRM_Name theCollisionWeightsName("collisionweights", "Collision Weights Field");
    static PRM_Name theCollisionVelocityName(GAS_NAME_COLLISIONVELOCITY, "Collision Velocity Field");
    static PRM_Name theSurfaceName(GAS_NAME_SURFACE, "Surface Field");
    static PRM_Default theSurfaceDefault(0, "surface");
    static PRM_Name theSurfaceWeightsName("surfaceweights", "Surface Weights Field");
    static PRM_Name theSurfacePressureName("surfacepressure", "Surface Pressure Field");
    static PRM_Default theSurfacePressureDefault(0, "surfacepressure");
    static PRM_Name thePressureName(GAS_NAME_PRESSURE, "Pressure Field");
    static PRM_Default thePressureDefault(0, "pressure");
    static PRM_Name theDensityName(GAS_NAME_DENSITY, "Density Field");
    static PRM_Name theMinDensityName("mindensity", "Min Density");
    static PRM_Name theMaxDensityName("maxdensity", "Max Density");
    static PRM_Default theMaxDensityDefault(100000);
    static PRM_Name theMinViscosityName("minviscosity", "Min Viscosity");
    static PRM_Default theMinViscosityDefault(0.01);
    static PRM_Name theValidName("valid", "Valid Field");
    static PRM_Name theScaleName(SIM_NAME_SCALE, "Scale");
    static PRM_Name theSupersamplingName("numsupersamples", "Samples Per Axis");
    static PRM_Name theFloatPrecisionName("floatprecision", "Float Precision");
    static PRM_Name theSchemeName("scheme", "Scheme");
    static PRM_Name theUseOpenCLName(SIM_NAME_OPENCL, "Use OpenCL");
    static PRM_Name theToleranceName(SIM_NAME_TOLERANCE, "Error Tolerance");
    static PRM_Default theToleranceDefault(1e-3);

    static PRM_Name theFloatPrecisionChoices[] = {
                                                PRM_Name("f32b", "Float 32 bit"),
                                                PRM_Name("f64b", "Float 64 bit"),
                                                PRM_Name(0)
                                            };

    static PRM_Name theSchemeChoices[] = {
                                        PRM_Name("stokes", "Full Stokes"),
                                        PRM_Name("decoupled", "Decoupled Viscosity -> Pressure"),
                                        PRM_Name("decoupledfancy", "Decoupled Pressure -> Viscosity -> Pressure"),
                                        PRM_Name("decouplednoexp", "Decoupled No-Expansion Viscosity -> Pressure"),
                                        PRM_Name("decouplednoexpfancy", "Decoupled No-Expansion Pressure -> Viscosity -> Pressure"),
                                        PRM_Name("pressureonly", "Pressure Only"),
                                        PRM_Name("viscosityonly", "Viscosity Only"),
                                        PRM_Name(0)
                                        };

    static PRM_ChoiceList theSchemeMenu(PRM_CHOICELIST_SINGLE, theSchemeChoices);
    static PRM_ChoiceList theFloatPrecisionMenu(PRM_CHOICELIST_SINGLE, theFloatPrecisionChoices);

    static PRM_Template    theTemplates[] = {
                                            PRM_Template(PRM_STRING, 1, &theVelocityName, &theVelocityDefault),
                                            PRM_Template(PRM_STRING, 1, &theViscosityName, &theViscosityDefault),
                                            PRM_Template(PRM_STRING, 1, &theSurfaceName, &theSurfaceDefault),
                                            PRM_Template(PRM_STRING, 1, &theSurfaceWeightsName),
                                            PRM_Template(PRM_STRING, 1, &theSurfacePressureName, &theSurfacePressureDefault),
                                            PRM_Template(PRM_STRING, 1, &theCollisionName, &theCollisionDefault),
                                            PRM_Template(PRM_STRING, 1, &theCollisionWeightsName),
                                            PRM_Template(PRM_STRING, 1, &theCollisionVelocityName),
                                            PRM_Template(PRM_STRING, 1, &thePressureName, &thePressureDefault),
                                            PRM_Template(PRM_STRING, 1, &theDensityName),
                                            PRM_Template(PRM_FLT, 1, &theMinDensityName, PRMoneDefaults, 0, 0, 0, &PRM_SpareData::unitsDensity),
                                            PRM_Template(PRM_FLT, 1, &theMaxDensityName, &theMaxDensityDefault, 0, 0, 0, &PRM_SpareData::unitsDensity),
                                            PRM_Template(PRM_FLT, 1, &theMinViscosityName, &theMinViscosityDefault),
                                            PRM_Template(PRM_STRING, 1, &theValidName),
                                            PRM_Template(PRM_FLT, 1, &theScaleName, PRMoneDefaults),
                                            PRM_Template(PRM_INT, 1, &theSupersamplingName, PRMtwoDefaults),
                                            PRM_Template(PRM_ORD, 1, &theFloatPrecisionName, PRMoneDefaults, &theFloatPrecisionMenu),
                                            PRM_Template(PRM_ORD, 1, &theSchemeName, PRMzeroDefaults, &theSchemeMenu),
                                            PRM_Template(PRM_TOGGLE, 1, &theUseOpenCLName, PRMoneDefaults),
                                            PRM_Template(PRM_FLT , 1, &theToleranceName, &theToleranceDefault),
                                            PRM_Template()
                                        };

    static SIM_DopDescription  theDopDescription(
                                                true,   // Should we make a DOP?
                                                "hdk_stokes",  // Internal name of the DOP.
                                                "Stokes",   // Label of the DOP
                                                "Solver",   // Default data name
                                                classname(),  // The type of this DOP, usually the class.
                                                theTemplates
                                                );  // Template list for generating the DOP
    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool SIM_Stokes::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) {


    return true;
}
