#include "SIM_Stokes.hpp"
#include "SolverCore.cpp"
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
#include <iostream>

using namespace Stokes;

void initializeSIM(void*) {
    Eigen::setNbThreads(12);
    #ifndef NDEBUG
    if ( Eigen::nbThreads() > 1 )
        std::cout << "Eigen is multithreaded" << std::endl;
    #endif
    IMPLEMENT_DATAFACTORY(SIM_Stokes);
}

SIM_Stokes::SIM_Stokes(const SIM_DataFactory* factory) : BaseClass(factory) {
}

SIM_Stokes::~SIM_Stokes() {
}

const SIM_DopDescription* SIM_Stokes::getDopDescription() {
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

    static PRM_Template theTemplates[] = {
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

    static SIM_DopDescription theDopDescription(
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



bool SIM_Stokes::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) {
    SIM_DataArray data;
    UT_StringArray datanames;

    /// =================== Import Required Fields ===================
    SIM_VectorField* velocity = getVectorField(obj, GAS_NAME_VELOCITY);
    const SIM_ScalarField* surface = getConstScalarField(obj, GAS_NAME_SURFACE);
    const SIM_ScalarField* collision = getConstScalarField(obj, GAS_NAME_COLLISION);


    /// =================== Check Required Fields ===================
    if (!velocity){
        addError(obj,SIM_MESSAGE, "No velocity detected", UT_ERROR_ABORT);
        return false;
    }
    if (!surface){
        addError(obj,SIM_MESSAGE, "No surface detected", UT_ERROR_ABORT);
        return false;
    }
    if (!collision){
        addError(obj,SIM_MESSAGE, "No collision surface detected", UT_ERROR_ABORT);
        return false;
    }
    if (!velocity->isFaceSampled()){
        addError(obj,SIM_MESSAGE, "Velocity field must be face sampled", UT_ERROR_ABORT);
        return false;
    }


    /// =================== Import Optional Fields ===================
    SIM_VectorField* valid = getVectorField(obj, "valid");
    const SIM_VectorField* collisionvel = getConstVectorField(obj, GAS_NAME_COLLISIONVELOCITY);
    const SIM_VectorField* colweights = getVectorField(obj, "collisionweights");
    const SIM_VectorField* surfweights = getVectorField(obj, "surfaceweights");
    const SIM_ScalarField* surfpressure = getConstScalarField(obj, "surfacepressure");
    //const SIM_ScalarField* pressure = getScalarField(obj, GAS_NAME_PRESSURE, true);
    const SIM_ScalarField* viscosity = getScalarField(obj, "viscosity");
    const SIM_ScalarField* density = getScalarField(obj, "density");


    /// =================== Check Optional Fields ===================
    if (!valid) {
        addError(obj,SIM_MESSAGE, "No valid field detected", UT_ERROR_MESSAGE);
    }
    if (valid && !valid->isAligned(velocity)){
        addError(obj,SIM_MESSAGE, "Valid field misaligned with velocity", UT_ERROR_ABORT);
        return false;
    }
    if (!surfpressure) {
        addError(obj,SIM_MESSAGE, "No surface pressure detected", UT_ERROR_MESSAGE);
    }
    if (!viscosity) {
        addError(obj,SIM_MESSAGE, "Viscosity field missing", UT_ERROR_WARNING);
    }
    if (!density) {
        addError(obj,SIM_MESSAGE, "Density field missing", UT_ERROR_WARNING);
    }


    /// =================== Get field configuration ===================
    fpreal dx = velocity->getVoxelSize(0).maxComponent();
    auto size = velocity->getSize();
    auto orig = velocity->getOrig();
    UT_Vector3 res = velocity->getTotalVoxelRes();
    exint nx = res.x(), ny = res.y(), nz = res.z();

    nx -= 1;
    ny -= 1;
    nz -= 1;
    std::cerr << " nx = " << nx << "; ny = " << ny << "; nz = " << nz << std::endl;


    /// =================== Validate Viscosity Field ===================
    SIM_RawField viscfielddata;
    SIM_RawField* viscfield = NULL;
    if ( viscosity ) {
        viscfield = const_cast<SIM_RawField*>(viscosity->getField());
    }else{
        viscfielddata.init(SIM_SAMPLE_CENTER,  orig, size, nx+1, ny+1, nz+1);
        viscfielddata.makeConstant(0);
        viscfield = & viscfielddata;
    }

    /// =================== Validate Density Field ===================
    SIM_RawField densfielddata;
    SIM_RawField* densfield = NULL;
    if ( density ) {
        densfield = const_cast<SIM_RawField*>(density->getField());
    }else{
        densfielddata.init(SIM_SAMPLE_CENTER,  orig, size, nx+1, ny+1, nz+1);
        densfielddata.makeConstant(1);
        densfield = & densfielddata;
    }

    assert( viscfield && densfield );

    fpreal scale = getScale();
    if ( SYSequalZero(scale) ) {
        return true; // no effect with zero scale
    }

    /// =================== Validate Collision Velocity Field ===================
    const SIM_RawField* colvel[3];
    SIM_RawField u_colvel, v_colvel, w_colvel;
    if (collisionvel) {
        colvel[0] = collisionvel->getField(0);
        colvel[1] = collisionvel->getField(1);
        colvel[2] = collisionvel->getField(2);
    }else{
        u_colvel.makeConstant(0);
        v_colvel.makeConstant(0);
        w_colvel.makeConstant(0);
        colvel[0] = & u_colvel;
        colvel[1] = & v_colvel;
        colvel[2] = & w_colvel;
    }


    /// =================== Validate Surface Pressure Field ===================
    const SIM_RawField* surfpres;
    SIM_RawField surfpresfield;
    if (surfpressure) {
        surfpres = surfpressure->getField();
    }else{
        surfpresfield.match(* surface->getField());
        surfpresfield.makeConstant(0);
        surfpres = & surfpresfield;
    }


    /// =================== Compute Volume Fraction Weights ===================
    SIM_RawField* surffield = const_cast<SIM_RawField*>(surface->getField());
    SIM_RawField* colfield = const_cast<SIM_RawField*>(collision->getField());

    SIM_RawField c_liquid_weights, u_liquid_weights, v_liquid_weights, w_liquid_weights;
    SIM_RawField xy_liquid_weights, xz_liquid_weights, yz_liquid_weights;
    SIM_RawField c_fluid_weights, u_fluid_weights, v_fluid_weights, w_fluid_weights;
    SIM_RawField xy_fluid_weights, xz_fluid_weights, yz_fluid_weights;

    // reuse face sampled weights if provided
    SIM_RawField* sweights[7] = {
        &c_liquid_weights,
        &xy_liquid_weights,
        &xz_liquid_weights,
        &yz_liquid_weights,
        0, 0, 0
    };

    SIM_RawField* cweights[7] = {
        &c_fluid_weights,
        &xy_fluid_weights,
        &xz_fluid_weights,
        &yz_fluid_weights,
        0, 0, 0
    };

    int ns = getNumSuperSamples();

    fpreal32 cval;
    bool is_surf_const = false;
    if ( surffield->field()->isConstant(&cval) && cval < 0) {
        is_surf_const = true;
    }

    bool is_col_const = false;
    if ( colfield->field()->isConstant(&cval) && cval < 0) {
        is_col_const = true;
    }

    { // Compute surfweights
        UT_PerfMonAutoSolveEvent event(this, "Compute Surface Weights");

        if ( surfweights ) {
            std::cerr << " surfweights->getField(0) = " << surfweights->getField(0) << "\n"
                      << " surfweights->getField(1) = " << surfweights->getField(1) << "\n"
                      << " surfweights->getField(2) = " << surfweights->getField(2) << "\n";
            sweights[4] = const_cast<SIM_RawField*>(surfweights->getField(0));
            sweights[5] = const_cast<SIM_RawField*>(surfweights->getField(1));
            sweights[6] = const_cast<SIM_RawField*>(surfweights->getField(2));
            for ( int i = 4; i < 7; ++i ){
                sweights[i]->setScaleDivideThreshold(1, NULL, NULL, MINWEIGHT);
            }
            std::cerr << " surfweights not eq null" << std::endl;
        }else{
            simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_FACEX,  ns, false, u_liquid_weights);
            simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_FACEY,  ns, false, v_liquid_weights);
            simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_FACEZ,  ns, false, w_liquid_weights);
            sweights[4] = &u_liquid_weights;
            sweights[5] = &v_liquid_weights;
            sweights[6] = &w_liquid_weights;
            std::cerr << " surfweights eq null" << std::endl;
        }

        simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_CENTER, ns, false, c_liquid_weights);
        simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_EDGEXY, ns, false, xy_liquid_weights);
        simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_EDGEXZ, ns, false, xz_liquid_weights);
        simEstimateVolumeFractions(surffield, is_surf_const, SIM_SAMPLE_EDGEYZ, ns, false, yz_liquid_weights);
    }

    { // Compute colweights
        UT_PerfMonAutoSolveEvent event(this, "Compute Collision Weights");

        if ( colweights ) {
            std::cerr << " colweights->getField(0) = " << colweights->getField(0) << "\n"
                      << " colweights->getField(1) = " << colweights->getField(1) << "\n"
                      << " colweights->getField(2) = " << colweights->getField(2) << "\n";
            cweights[4] = const_cast<SIM_RawField*>(colweights->getField(0));
            cweights[5] = const_cast<SIM_RawField*>(colweights->getField(1));
            cweights[6] = const_cast<SIM_RawField*>(colweights->getField(2));
            for ( int i = 4; i < 7; ++i ){
                cweights[i]->setScaleDivideThreshold(1, NULL, NULL, MINWEIGHT);
            }
            std::cerr << " colweights not eq null" << std::endl;
        }else{
            simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_FACEX,  ns, false, u_fluid_weights);
            simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_FACEY,  ns, false, v_fluid_weights);
            simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_FACEZ,  ns, false, w_fluid_weights);
            cweights[4] = &u_fluid_weights;
            cweights[5] = &v_fluid_weights;
            cweights[6] = &w_fluid_weights;
            std::cerr << " colweights eq null" << std::endl;
        }

        simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_CENTER, ns, false, c_fluid_weights);
        simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_EDGEXY, ns, false, xy_fluid_weights);
        simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_EDGEXZ, ns, false, xz_fluid_weights);
        simEstimateVolumeFractions(colfield, is_col_const, SIM_SAMPLE_EDGEYZ, ns, false, yz_fluid_weights);
    }

    //#ifndef NDEBUG
        for (int i = 0; i < 7; ++i ) {
            if (!(sweights[i] && cweights[i])) {

                addError(obj,SIM_MESSAGE, "Did not match sweights value between cweights value.", UT_ERROR_ABORT);
                return false;
            }
        } // make sure we got all of them

    //#endif
    /// =================== End of Computing Volume Fraction Weights ===================


    FloatPrecision float_precision( getFloatPrecision() );
    SolverResult result = NOCHANGE;


    /// =================== Solve System and update Velocities ===================
    if ( float_precision == FLOAT32 ) {
        sim_stokesSolver<fpreal32> solver(*this, obj, nx, ny, nz, dx, timestep);
        solver.classifyAndBuildIndices(sweights, cweights);
        result = solver.solve(*surffield, sweights, cweights, *viscfield, *densfield, colvel, *surfpres, valid, *velocity);
    }else{
        assert( float_precision == FLOAT64 ); // only one option left
        sim_stokesSolver<fpreal64> solver(*this, obj, nx, ny, nz, dx, timestep);
        solver.classifyAndBuildIndices(sweights, cweights);
        result = solver.solve(*surffield, sweights, cweights, *viscfield, *densfield, colvel, *surfpres, valid, *velocity);
    }

    if ( result == SUCCESS ) {
        velocity->pubHandleModification();
        if ( valid ){
            valid->pubHandleModification();
        }
    }

    return result == SUCCESS || result == NOCHANGE;
}
