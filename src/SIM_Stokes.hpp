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
                GET_DATA_FUNC_F(SIM_NAME_SCALE, Scale);
                GET_DATA_FUNC_E("scheme", Scheme, Stokes::Scheme);
                GET_DATA_FUNC_I("numsupersamples", NumSuperSamples);
                GET_DATA_FUNC_E("floatprecision", FloatPrecision, Stokes::FloatPrecision);
                GET_DATA_FUNC_B(SIM_NAME_OPENCL, UseOpenCL);
                GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, Tolerance);
                GET_DATA_FUNC_F("maxdensity", MaxDensity);
                GET_DATA_FUNC_F("mindensity", MinDensity);
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

                DECLARE_STANDARD_GETCASTTOTYPE();
                DECLARE_DATAFACTORY(SIM_Stokes,
                                    GAS_SubSolver,
                                    "Stokes",
                                    getDopDescription()
                                    );
    };

} // End namespace

#endif
