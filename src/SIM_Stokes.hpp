#ifndef __SIM_Stokes_hpp__
#define __SIM_Stokes_hpp__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

namespace Stokes {

    enum Scheme {
        STOKES,
        DECOUPLED,
        DECOUPLED_FANCY,
        DECOUPLED_NOEXPANSION,
        DECOUPLED_NOEXPANSION_FANCY,
        PRESSURE_ONLY,
        VISCOSITY_ONLY
    };

    enum FloatPrecision {
        FLOAT32,
        FLOAT64
    };

    class SIM_Stokes : public GAS_SubSolver {
        public:
                //define `SIM_Stokes::getScale() -> const fpreal64`
                GET_DATA_FUNC_F(SIM_NAME_SCALE, Scale);

                //define `SIM_Stokes::getScheme() -> const Stokes::Scheme`
                GET_DATA_FUNC_E("scheme", Scheme, Stokes::Scheme);

                //define `SIM_Stokes::getNumSuperSamples() -> const int64`
                GET_DATA_FUNC_I("numsupersamples", NumSuperSamples);

                //define `SIM_Stokes::getFloatPrecision() -> const Stokes::FloatPrecision`
                GET_DATA_FUNC_E("floatprecision", FloatPrecision, Stokes::FloatPrecision);

                //define `SIM_Stokes::getUseOpenCL() -> const bool`
                GET_DATA_FUNC_B(SIM_NAME_OPENCL, UseOpenCL);

                //define `SIM_Stokes::getTolerance() -> const fpreal64`
                GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, Tolerance);

                //define `SIM_Stokes::getMaxDensity() -> const fpreal64`
                GET_DATA_FUNC_F("maxdensity", MaxDensity);

                //define `SIM_Stokes::getMinDensity() -> const fpreal64`
                GET_DATA_FUNC_F("mindensity", MinDensity);

                //define `SIM_Stokes::getMinViscosity() -> const fpreal64`
                GET_DATA_FUNC_F("minviscosity", MinViscosity);


            protected:
                explicit  SIM_Stokes(const SIM_DataFactory *factory);
                virtual  ~SIM_Stokes();

                bool shouldMultiThread(const SIM_RawField *field) const {
                    return field->field()->numTiles() > 1;
                }

                bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) override;

            private:
                static const SIM_DopDescription* getDopDescription();

                // define `protected: SIM_Stokes::getCastToType(const UT_StringRef&) -> const override void*`
                DECLARE_STANDARD_GETCASTTOTYPE();

                // define `SIM_StokesFactory : SIM_DataFactory` class
                DECLARE_DATAFACTORY(SIM_Stokes,
                                    GAS_SubSolver,
                                    "Stokes",
                                    getDopDescription()
                                    );
    };

} // End namespace

#endif
