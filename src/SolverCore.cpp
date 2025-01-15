#include "SIM_Stokes.hpp"
#include "UT/UT_SparseMatrix.h"
#include "util/eigen.h"
#include <UT/UT_DSOVersion.h>
#include <UT/UT_PerfMonAutoEvent.h>
#include <UT/UT_VoxelArray.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_FieldSampler.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_RawIndexField.h>
#include <SIM/SIM_Object.h>
#include <CE/CE_Vector.h>
#include <CE/CE_SparseMatrix.h>

namespace Stokes{
    enum FieldIndex {
        CENTER = 0,
        EDGEXY = 1,
        EDGEXZ = 2,
        EDGEYZ = 3,
        FACEX = 4,
        FACEY = 5,
        FACEZ = 6
    };
    enum SolveType {
        COLLISION = -3,
        AIR = -2,
        INVALIDIDX = -1,
        SOLVED = 1,
    };
    enum SolverResult{
        NOCONVERGE = 0,
        SUCCESS = 1,
        NOCHANGE = 2,
        FAILED = -1,
        INVALID = -2
    };


    // import all enum values for convenience
    using namespace Stokes;


    template<typename T> class sim_stokesSolver {
        using MatrixType = UT_SparseMatrixELLT<T, /*colmajor*/true>;
        using VectorType = UT_VectorT<T>;
        using BlockMatrixType = Eigen::SparseMatrix<T>;
        using BlockVectorType = VecX<T>;

        public:

            // Structs
            struct sim_buildSystemParms {
                const UT_VoxelArrayF &c_vol_liquid;
                const UT_VoxelArrayF &ez_vol_liquid;
                const UT_VoxelArrayF &ey_vol_liquid;
                const UT_VoxelArrayF &ex_vol_liquid;
                const UT_VoxelArrayF &u_vol_liquid;
                const UT_VoxelArrayF &v_vol_liquid;
                const UT_VoxelArrayF &w_vol_liquid;

                const UT_VoxelArrayF &c_vol_fluid;
                const UT_VoxelArrayF &ez_vol_fluid;
                const UT_VoxelArrayF &ey_vol_fluid;
                const UT_VoxelArrayF &ex_vol_fluid;
                const UT_VoxelArrayF &u_vol_fluid;
                const UT_VoxelArrayF &v_vol_fluid;
                const UT_VoxelArrayF &w_vol_fluid;

                const UT_VoxelArrayF &u;
                const UT_VoxelArrayF &v;
                const UT_VoxelArrayF &w;

                const UT_VoxelArrayF &u_solid;
                const UT_VoxelArrayF &v_solid;
                const UT_VoxelArrayF &w_solid;

                const UT_VoxelArrayF &viscosity;
                const UT_VoxelArrayF &density;
                const UT_VoxelArrayF &surfpres;

                fpreal minrho;
                fpreal maxrho;
            };

            struct sim_updateVelocityParms {
                sim_updateVelocityParms(const SIM_RawField * const* sweights,
                                        const SIM_RawField * const* solid_vel,
                                        const SIM_RawField & densfield,
                                        const SIM_RawField & surfpres,
                                        fpreal min_density,
                                        fpreal max_density
                                        )
                : c_vol_liquid( *sweights[0]->field())
                , ez_vol_liquid(*sweights[1]->field())
                , ey_vol_liquid(*sweights[2]->field())
                , ex_vol_liquid(*sweights[3]->field())
                , u_vol_liquid( *sweights[4]->field())
                , v_vol_liquid( *sweights[5]->field())
                , w_vol_liquid( *sweights[6]->field())

                , u_solid(*solid_vel[0]->field())
                , v_solid(*solid_vel[1]->field())
                , w_solid(*solid_vel[2]->field())

                , density(*densfield.field())
                , surfpres(*surfpres.field())

                , minrho(min_density)
                , maxrho(max_density)
              {
              }
                const UT_VoxelArrayF &c_vol_liquid;
                const UT_VoxelArrayF &ez_vol_liquid;
                const UT_VoxelArrayF &ey_vol_liquid;
                const UT_VoxelArrayF &ex_vol_liquid;
                const UT_VoxelArrayF &u_vol_liquid;
                const UT_VoxelArrayF &v_vol_liquid;
                const UT_VoxelArrayF &w_vol_liquid;

                const UT_VoxelArrayF &u_solid;
                const UT_VoxelArrayF &v_solid;
                const UT_VoxelArrayF &w_solid;

                const UT_VoxelArrayF &density;
                const UT_VoxelArrayF &surfpres;

                fpreal minrho;
                fpreal maxrho;
            };

            // Constructor
            sim_stokesSolver(SIM_Stokes& solver, SIM_Object* obj, int nx, int ny, int nz, float dx, float dt)  : ni(nx), nj(ny), nk(nz), dx(dx), dt(dt)
                , myNumStressVars(0)
                , myNumVelocityVars(0)
                , myNumPressureVars(0)
                , myCollisionIndex(std::numeric_limits<int>::max())
                , mySolver(solver)
                , myObject(obj) {
                // pass
            }

            //***********************************************************************************//
            /// Main entry point into the solver. This builds the system, solves it and
            /// updates velocities
            //***********************************************************************************//
            SolverResult solve( const SIM_RawField & phi,
                                const SIM_RawField * const* sweights,
                                const SIM_RawField * const* cweights,
                                const SIM_RawField & viscosity,
                                const SIM_RawField & density,
                                const SIM_RawField * const* solid_vel,
                                const SIM_RawField & surfpres,
                                SIM_VectorField * valid,
                                SIM_VectorField & vel
                                ) const;

            SolverResult solveStokes(   const SIM_RawField * const* sweights,
                                        const SIM_RawField * const* cweights,
                                        const SIM_RawField & viscosity,
                                        const SIM_RawField & density,
                                        const SIM_RawField * const* solid_vel,
                                        const SIM_RawField & surfpres,
                                        SIM_VectorField * valid,
                                        SIM_VectorField & vel
                                        ) const;

            void buildSystem(   MatrixType &A,
                                VectorType &b,
                                const sim_buildSystemParms& parms
                                ) const;

            // helper for solveStokes
            SolverResult solveSystem(   const MatrixType &A,
                                        const VectorType &b,
                                        VectorType &x,
                                        bool use_opencl
                                        ) const;

            bool isInSystem(exint idx) const {
                return idx >= 0;
            }

            // return true if velocity index represents a collision velocity in the system
            bool isCollision(exint idx) const {
                return idx == COLLISION || idx >= myCollisionIndex;
            }

            template<int AXIS> auto ghostFluidSurfaceTensionPressure(int i, int j, int k, float uweight, const UT_VoxelArrayF & sp) const -> T;


            /// =================== Main Func and Multithreded ===================
            void classifyAndBuildIndices(const SIM_RawField * const* surf_weights, const SIM_RawField * const* col_weights);

            void initAndClassifyIndex(   const SIM_RawField * const* surf_weights,
                                            const SIM_RawField * const* col_weights,
                                            SIM_RawIndexField &index,
                                            FieldIndex fidx
                                        );

            void buildIndex(SIM_RawIndexField &index,
                                FieldIndex fidx,
                                exint &maxindex
                            );

            void buildVelocityIndices(const SIM_RawField * const* surf_weights, const SIM_RawField * const* col_weights);

            void buildCollisionIndex(   SIM_RawIndexField &index,
                                        FieldIndex fidx,
                                        exint &maxindex
                                        );

            SolveType solveType(const SIM_RawField * const* surf_weights,
                                const SIM_RawField * const* col_weights,
                                int i, int j, int k,
                                FieldIndex fidx
                                ) const;

            // build member index fields
            THREADED_METHOD4(sim_stokesSolver, index.shouldMultiThread(),
                           classifyIndexField,
                           const SIM_RawField * const*, surf_weights,
                           const SIM_RawField * const*, col_weights,
                           SIM_RawIndexField &, index,
                           FieldIndex, fidx);

            void classifyIndexFieldPartial( const SIM_RawField * const* surf_weights,
                                            const SIM_RawField * const* col_weights,
                                            SIM_RawIndexField &index,
                                            FieldIndex fidx,
                                            const UT_JobInfo &info
                                        );


            /// =================== Velocity Update func and Multithreded ===================
            THREADED_METHOD5_CONST( sim_stokesSolver, vel.getField(axis)->shouldMultiThread(),
                                    updateVelocities,
                                    const VectorType &, x,
                                    const sim_updateVelocityParms &, parms,
                                    SIM_VectorField *, valid, // output is explicit
                                    SIM_VectorField &, vel,
                                    int, axis
                                    )

             void updateVelocitiesPartial(  const VectorType &x,
                                            const sim_updateVelocityParms &parms,
                                            SIM_VectorField *valid,
                                            SIM_VectorField &vel,
                                            int axis,
                                            const UT_JobInfo &info
                                            ) const;

            /// =================== Other Func and Multithreded ===================
            THREADED_METHOD3_CONST( sim_stokesSolver, myCentralIndex.shouldMultiThread(),
                                    addCenterTerms,
                                    MatrixType&, A,
                                    VectorType&, b,
                                    const sim_buildSystemParms&, parms
                                    );
            void addCenterTermsPartial( MatrixType& A,
                                        VectorType &b,
                                        const sim_buildSystemParms& parms,
                                        const UT_JobInfo& info
                                        ) const;
            THREADED_METHOD3_CONST( sim_stokesSolver, myTxyIndex.shouldMultiThread(),
                                    addTxyTerms,
                                    MatrixType&, A,
                                    VectorType&, b,
                                    const sim_buildSystemParms&, parms
                                    );
            void addTxyTermsPartial(MatrixType& A,
                                    VectorType &b,
                                    const sim_buildSystemParms& parms,
                                    const UT_JobInfo& info
                                    ) const;
            THREADED_METHOD3_CONST( sim_stokesSolver, myTxzIndex.shouldMultiThread(),
                                    addTxzTerms,
                                    MatrixType&, A,
                                    VectorType&, b,
                                    const sim_buildSystemParms&, parms
                                    );
            void addTxzTermsPartial(MatrixType& A,
                                    VectorType &b,
                                    const sim_buildSystemParms& parms,
                                    const UT_JobInfo& info
                                    ) const;
            THREADED_METHOD3_CONST( sim_stokesSolver, myTyzIndex.shouldMultiThread(),
                                    addTyzTerms,
                                    MatrixType&, A,
                                    VectorType&, b,
                                    const sim_buildSystemParms&, parms
                                    );
            void addTyzTermsPartial(MatrixType& A,
                                    VectorType &b,
                                    const sim_buildSystemParms& parms,
                                    const UT_JobInfo& info
                                    ) const;


        private:
            // date member
            int                 ni, nj, nk;
            float               dx, dt;
            int                 myNumPressureVars;
            int                 myNumVelocityVars; // including collision vars
            int                 myNumStressVars;
            int                 myCollisionIndex; // first velocity collision index (used in decoupled and blockwise solves)
            SIM_Stokes&         mySolver;
            SIM_Object*         myObject;
            SIM_RawIndexField   myCentralIndex, myTxyIndex, myTxzIndex, myTyzIndex; // stokes system indices

            // additional index fields for decoupled systems (for Stokes these just act as
            // SolveType flags since velocities don't get an actual index)
            SIM_RawIndexField myUIndex, myVIndex, myWIndex;


            // out of bounds checks
            bool c_oob(int i, int j, int k) const {
                return i < 0 || i > ni-1 || j < 0 || j > nj-1 || k < 0 || k > nk-1;
            }
            bool tyz_oob(int i, int j, int k) const {
                return i < 0 || i > ni-1 || j < 0 || j > nj || k < 0 || k > nk;
            }
            bool txz_oob(int i, int j, int k) const {
                return i < 0 || i > ni || j < 0 || j > nj-1 || k < 0 || k > nk;
            }
            bool txy_oob(int i, int j, int k) const {
                return i < 0 || i > ni || j < 0 || j > nj || k < 0 || k > nk-1;
            }
            bool u_oob(int i, int j, int k) const {
                return i < 0 || i > ni || j < 0 || j > nj-1 || k < 0 || k > nk-1;
            }
            bool v_oob(int i, int j, int k) const {
                return i < 0 || i > ni-1 || j < 0 || j > nj || k < 0 || k > nk-1;
            }
            bool w_oob(int i, int j, int k) const {
                return i < 0 || i > ni-1 || j < 0 || j > nj-1 || k < 0 || k > nk;
            }

            // System index getters
            // NOTE: to save on indirection, we use c_index to store the first 3 indices:
            // 1 for pressure, and 2 for txx and tyy respectively. Thus for instance the
            // index of tyy is myCentralIndex(i,j,k) + 2*solver.myNumPressureVars. As a result
            // NOTE: we do an out of bounds check because it should be faster than doing the more
            // general .getValue() call
            exint p_idx(int i, int j, int k) const {
                return c_oob(i,j,k) ? exint(INVALIDIDX) : myCentralIndex(i,j,k);
            }

            exint txx_idx(int i, int j, int k) const {
                return c_oob(i,j,k) ? exint(INVALIDIDX) : (myCentralIndex(i,j,k) + (isInSystem(myCentralIndex(i,j,k)) ? myNumPressureVars : 0));
            }
            exint tyy_idx(int i, int j, int k) const {
                return c_oob(i,j,k) ? exint(INVALIDIDX) : (myCentralIndex(i,j,k) + (isInSystem(myCentralIndex(i,j,k)) ? 2*myNumPressureVars : 0));
            }
            exint tyz_idx(int i, int j, int k) const {
                return tyz_oob(i,j,k) ? exint(INVALIDIDX) : myTyzIndex(i,j,k);
            }
            exint txz_idx(int i, int j, int k) const {
                return txz_oob(i,j,k) ? exint(INVALIDIDX) : myTxzIndex(i,j,k);
            }
            exint txy_idx(int i, int j, int k) const {
                return txy_oob(i,j,k) ? exint(INVALIDIDX) : myTxyIndex(i,j,k);
            }

    }; // End of sim_stokesSolver Class

} // End of namespace

using namespace Stokes;
/// =================== Scale Field Multithreded ===================
struct FieldArithmetic {
    THREADED_METHOD2_CONST( FieldArithmetic, A.shouldMultiThread(),
                            scale,
                            SIM_RawField&, A,
                            fpreal, scale
                            )

    void scalePartial(  SIM_RawField& A,
                        fpreal scale,
                        const UT_JobInfo& info
                        ) const;
};

void FieldArithmetic::scalePartial(SIM_RawField& A, fpreal scale, const UT_JobInfo &info) const {
    // compute A = scale*A;
    UT_VoxelArrayIteratorF vit;
    A.getPartialRange(vit, info);
    vit.setCompressOnExit(true);
    vit.detectInterrupts();
    auto op = [&scale](fpreal32 a) { return a * scale; };
    vit.applyOperation(op);
}


/// =================== Implement sim_stokesSolver::solve() ===================
// Note: Valid field is optional. It specifies which velocity samples were updated
template<typename T> SolverResult sim_stokesSolver<T>::solve(   const SIM_RawField & surf,
                                                                const SIM_RawField * const* sweights,
                                                                const SIM_RawField * const* cweights,
                                                                const SIM_RawField & viscosity,
                                                                const SIM_RawField & density,
                                                                const SIM_RawField * const* colvel,
                                                                const SIM_RawField & surfpres,
                                                                SIM_VectorField * valid,
                                                                SIM_VectorField & vel
                                                                ) const {
    SolverResult result = NOCHANGE;
    result = solveStokes(sweights, cweights, viscosity, density, colvel, surfpres, valid, vel);
    return result;
}

/// =================== Implement sim_stokesSolver::solveStokes() ===================
// Chain `solve() -> solveStokes()`
template<typename T> SolverResult sim_stokesSolver<T>::solveStokes( const SIM_RawField * const* sweights,
                                                                    const SIM_RawField * const* cweights,
                                                                    const SIM_RawField & viscfield,
                                                                    const SIM_RawField & densfield,
                                                                    const SIM_RawField * const* solid_vel,
                                                                    const SIM_RawField & surfpres,
                                                                    SIM_VectorField * valid,
                                                                    SIM_VectorField & vel
                                                                    ) const {
    sim_buildSystemParms parms{ *sweights[0]->field(),
                                *sweights[1]->field(),
                                *sweights[2]->field(),
                                *sweights[3]->field(),
                                *sweights[4]->field(),
                                *sweights[5]->field(),
                                *sweights[6]->field(),

                                *cweights[0]->field(),
                                *cweights[1]->field(),
                                *cweights[2]->field(),
                                *cweights[3]->field(),
                                *cweights[4]->field(),
                                *cweights[5]->field(),
                                *cweights[6]->field(),

                                *vel.getField(0)->field(),
                                *vel.getField(1)->field(),
                                *vel.getField(2)->field(),

                                *solid_vel[0]->field(),
                                *solid_vel[1]->field(),
                                *solid_vel[2]->field(),

                                *viscfield.field(),
                                *densfield.field(),
                                *surfpres.field(),

                                mySolver.getMinDensity(),
                                mySolver.getMaxDensity()
    };

    auto system_size = myNumPressureVars + myNumStressVars;
    // 29 is the max non zeros per row in the stokes system
    MatrixType A(system_size, 29);
    VectorType b(0, system_size-1);
    VectorType x(0, system_size-1);

    // Build the Main Stokes System
    buildSystem(A, b, parms);

    #ifndef NDEBUG
        for (int row = 0; row < A.getNumRows(); ++row) {
            auto idx = A.index(row,0);
            auto val = A.getColumns()[idx];
            assert( val != -1 );
        }
    #endif

    auto result = solveSystem(A, b, x, mySolver.getUseOpenCL());
    if (result == SUCCESS) {
        UT_PerfMonAutoSolveEvent event(&mySolver, "Update Velocity");
        sim_updateVelocityParms parms(  sweights, solid_vel, densfield, surfpres,
                                        mySolver.getMinDensity(), mySolver.getMaxDensity()
                                        );
        for ( int axis = 0; axis < 3; ++axis ) {
            if ( valid ) {
                valid->getField(axis)->makeConstant(0);
            }
            updateVelocities(x, parms, valid, vel, axis);
        }
  }
  return result;
}


/// =================== Implement sim_stokesSolver::buildSystem() ===================
// Chain `solve() -> solveStokes() -> buildSystem()`
template<typename T> void sim_stokesSolver<T>::buildSystem( MatrixType &A,
                                                            VectorType &b,
                                                            const sim_buildSystemParms& parms
                                                            ) const {
    UT_PerfMonAutoSolveEvent event(&mySolver, "Build System");
    b.zero();
    addCenterTerms(A,b,parms);
    addTxyTerms(A,b,parms);
    addTxzTerms(A,b,parms);
    addTyzTerms(A,b,parms);
    //A.sortRows(); // not necessary since we do this manually?
}


/// =================== Implement sim_stokesSolver::solveSystem() ===================
// Chain `solve() -> solveStokes() -> solveSystem()`
template<typename T> SolverResult sim_stokesSolver<T>::solveSystem( const MatrixType &A,
                                                                    const VectorType &b,
                                                                    VectorType &x,
                                                                    bool use_opencl
                                                                    ) const {
    b.testForNan();
    auto system_size = b.length();
    assert( b.length() == A.getNumRows() );
    if ( !system_size ) {
        return NOCHANGE;
    }

    UT_PerfMonAutoSolveEvent event(&mySolver, "Solve Stokes Linear System");

    T tol = mySolver.getTolerance();

    int iterations = 0; // report these later
    float error = 0;
    #ifndef CE_ENABLED
        if (use_opencl) {
            mySolver.addError(myObject, SIM_NO_OPENCL, 0, UT_ERROR_ABORT);
            return FAILED;
        }
    #else
        if (use_opencl) {
            CE_Context *context = CE_Context::getContext();
            use_opencl = !context->isCPU();
        }
        x.zero();
        if (use_opencl) {
            try {
                GAS_ScopedOCLErrorSink  errorsink(myObject, &mySolver);
                CE_SparseMatrixELLT<T>  Ac;
                CE_VectorT<T>           xc, bc;
                Ac.initFromMatrix(A);
                xc.initFromVector(x);
                bc.initFromVector(b);
                error = Ac.solveConjugateGradient(xc, bc, tol, system_size*3, &iterations);
                xc.matchAndCopyToVector(x);
            } catch (cl::Error &err) {
                mySolver.addError(myObject, SIM_MESSAGE, "No velocity detected", UT_ERROR_ABORT);
                return FAILED;
            }
        } else {

        }
    #endif
    error = A.solveConjugateGradient(x, b, NULL, tol, system_size*3, &iterations);

    UT_WorkBuffer extra_info;
    extra_info.sprintf("Iterations=%d, Error=%.6f", int(iterations), error);
    event.setExtraInfo(extra_info.buffer());
    return SUCCESS;
}


/// =================== Implement sim_stokesSolver::updateVelocitiesPartial() ===================
// Chain `solve() -> solveStokes() -> updateVelocities() -> updateVelocitiesPartial()`
template<typename T> void sim_stokesSolver<T>::updateVelocitiesPartial( const VectorType &x,
                                                                        const sim_updateVelocityParms & parms,
                                                                        SIM_VectorField * valid,
                                                                        SIM_VectorField &vel,
                                                                        int axis,
                                                                        const UT_JobInfo &info
                                                                        ) const {
    // Update velocities based on the pressures and stresses determined by the solver.
    UT_VoxelArrayF &u = *vel.getField(axis)->fieldNC();

    // edge-centred quantities
    auto txy = [&](int i, int j, int k) { return !isInSystem(txy_idx(i,j,k)) ? 0 : x(txy_idx(i,j,k)); };
    auto txz = [&](int i, int j, int k) { return !isInSystem(txz_idx(i,j,k)) ? 0 : x(txz_idx(i,j,k)); };
    auto tyz = [&](int i, int j, int k) { return !isInSystem(tyz_idx(i,j,k)) ? 0 : x(tyz_idx(i,j,k)); };

    // cell centered quantities
    auto txx = [&](int i, int j, int k) { return !isInSystem(txx_idx(i,j,k)) ? 0 : x(txx_idx(i,j,k)); };
    auto tyy = [&](int i, int j, int k) { return !isInSystem(tyy_idx(i,j,k)) ? 0 : x(tyy_idx(i,j,k)); };
    auto p   = [&](int i, int j, int k) { return !isInSystem(p_idx(i,j,k)  ) ? 0 : x(p_idx(i,j,k)); };

    if ( axis == 0 ) {
        UT_VoxelProbeAverage<float,-1,0,0> rhox;
        rhox.setArray(&parms.density);
        UT_VoxelArrayIteratorF vit(&u);
        vit.splitByTile(info);

        for ( vit.rewind(); !vit.atEnd(); vit.advance() ) {
            int i = vit.x(), j = vit.y(), k = vit.z();
            int idx = myUIndex(i,j,k);
            if (isCollision(idx)) {
                vit.setValue(parms.u_solid.getValue(i,j,k));
                if (valid) {
                    valid->getField(axis)->fieldNC()->setValue(i,j,k,1);
                }
            } else if (isInSystem(idx)) {
                if (valid) {
                    valid->getField(axis)->fieldNC()->setValue(i,j,k,1);
                    auto gfp = ghostFluidSurfaceTensionPressure<0>(i,j,k, parms.u_vol_liquid(i,j,k), parms.surfpres);
                    rhox.setIndex(vit);
                    auto rho = SYSclamp(rhox.getValue(), parms.minrho, parms.maxrho);
                    auto factor = dt / (dx * rho * parms.u_vol_liquid(i,j,k));
                    // pressure
                    vit.setValue(u(i,j,k) + factor * (parms.c_vol_liquid.getValue(i-1,j,k)*p(i-1,j,k) - parms.c_vol_liquid.getValue(i,j,k)*p(i,j,k)
                        // stress
                        + ((parms.c_vol_liquid.getValue(i,j,k)    *txx(i,j,k)   - parms.c_vol_liquid.getValue(i-1,j,k) *txx(i-1,j,k))
                        +  (parms.ez_vol_liquid.getValue(i,j+1,k) *txy(i,j+1,k) - parms.ez_vol_liquid.getValue(i,j,k)  *txy(i,j,k))
                        +  (parms.ey_vol_liquid.getValue(i,j,k+1) *txz(i,j,k+1) - parms.ey_vol_liquid.getValue(i,j,k)  *txz(i,j,k))))
                        - factor * gfp
                        );
                }
            } else {
                vit.setValue(0);
            }
        }
    } else if ( axis == 1 ) {
        UT_VoxelProbeAverage<float,0,-1,0> rhoy;
        rhoy.setArray(&parms.density);
        UT_VoxelArrayIteratorF vit(&u);
        vit.splitByTile(info);

        for ( vit.rewind(); !vit.atEnd(); vit.advance() ) {
            int i = vit.x(), j = vit.y(), k = vit.z();
            int idx = myVIndex(i,j,k);
            if (isCollision(idx)) {
                vit.setValue(parms.v_solid.getValue(i,j,k));
                if (valid)
                valid->getField(axis)->fieldNC()->setValue(i,j,k,1);
            } else if (isInSystem(idx)) {
                if (valid) {
                valid->getField(axis)->fieldNC()->setValue(i,j,k,1);
                auto gfp = ghostFluidSurfaceTensionPressure<1>(i,j,k, parms.v_vol_liquid(i,j,k), parms.surfpres);
                rhoy.setIndex(vit);
                auto rho = SYSclamp(rhoy.getValue(), parms.minrho, parms.maxrho);
                auto factor = dt / (dx * rho * parms.v_vol_liquid(i,j,k));
                //pressure
                vit.setValue(u(i,j,k) + factor * (parms.c_vol_liquid.getValue(i,j-1,k)*p(i,j-1,k) - parms.c_vol_liquid.getValue(i,j,k)*p(i,j,k)
                    //stress
                    + ((parms.ez_vol_liquid.getValue(i+1,j,k) *txy(i+1,j,k) - parms.ez_vol_liquid.getValue(i,j,k)  *txy(i,j,k))
                    +  (parms.c_vol_liquid.getValue(i,j,k)    *tyy(i,j,k)   - parms.c_vol_liquid.getValue(i,j-1,k) *tyy(i,j-1,k))
                    +  (parms.ex_vol_liquid.getValue(i,j,k+1) *tyz(i,j,k+1) - parms.ex_vol_liquid.getValue(i,j,k)  *tyz(i,j,k))))
                    - factor * gfp
                    );
                }
            } else {
                vit.setValue(0);
            }
        }
    } else if ( axis == 2 ) {
        UT_VoxelProbeAverage<float,0,0,-1> rhoz;
        rhoz.setArray(&parms.density);
        UT_VoxelArrayIteratorF vit(&u);
        vit.splitByTile(info);

        for ( vit.rewind(); !vit.atEnd(); vit.advance() ) {
            int i = vit.x(), j = vit.y(), k = vit.z();
            int idx = myWIndex(i,j,k);
            if (isCollision(idx)) {
                vit.setValue(parms.w_solid.getValue(i,j,k));
                if (valid){
                    valid->getField(axis)->fieldNC()->setValue(i,j,k,1);
                }
            } else if (isInSystem(idx)) {
                if (valid) {
                valid->getField(axis)->fieldNC()->setValue(i,j,k,1);

                auto gfp = ghostFluidSurfaceTensionPressure<2>(i,j,k, parms.w_vol_liquid(i,j,k), parms.surfpres);
                rhoz.setIndex(vit);
                auto rho = SYSclamp(rhoz.getValue(), parms.minrho, parms.maxrho);
                auto factor =  dt / (dx * rho * parms.w_vol_liquid(i,j,k));
                //pressure
                vit.setValue(u(i,j,k) + factor * (parms.c_vol_liquid.getValue(i,j,k-1)*p(i,j,k-1) - parms.c_vol_liquid.getValue(i,j,k)*p(i,j,k)
                    //stress
                    + ((parms.ey_vol_liquid.getValue(i+1,j,k)*txz(i+1,j,k) - parms.ey_vol_liquid.getValue(i,j,k)  *txz(i,j,k))
                    +  (parms.ex_vol_liquid.getValue(i,j+1,k)*tyz(i,j+1,k) - parms.ex_vol_liquid.getValue(i,j,k)  *tyz(i,j,k))
                    -  (parms.c_vol_liquid.getValue(i,j,k)   *txx(i,j,k)   - parms.c_vol_liquid.getValue(i,j,k-1) *txx(i,j,k-1))
                    -  (parms.c_vol_liquid.getValue(i,j,k)   *tyy(i,j,k)   - parms.c_vol_liquid.getValue(i,j,k-1) *tyy(i,j,k-1))))
                    - factor * gfp
                    );
                }
            } else {
                vit.setValue(0);
            }
        }
    } else{
        assert(0); // uknown dimension
    }
}


/// =================== Implement sim_stokesSolver::ghostFluidSurfaceTensionPressure() ===================
// Chain `solve() -> solveStokes() -> updateVelocities() -> updateVelocitiesPartial() -> ghostFluidSurfaceTensionPressure()`
// POST: if there the velocity sample lies at the boundary, return the
// appropriate ghost fluid pressure at the air-liquid interface, which lies
// within the given velocity voxel
template<typename T> template<int AXIS> auto sim_stokesSolver<T>::ghostFluidSurfaceTensionPressure( int i, int j, int k, float uweight,
                                                                                                    const UT_VoxelArrayF & sp
                                                                                                    ) const -> T {
    if ( !uweight ) return 0;
    exint uidx = -1;
    exint pidx0 = -1;
    exint pidx1 = myCentralIndex(i,j,k);
    auto p1 = sp.getValue(i,j,k);
    float p0 = 0;
    switch (AXIS) {
        case 0:
            uidx = myUIndex(i,j,k);
            pidx0 = p_idx(i-1,j,k);
            p0 = sp.getValue(i-1,j,k);
            break;
        case 1:
            uidx = myVIndex(i,j,k);
            pidx0 = p_idx(i,j-1,k);
            p0 = sp.getValue(i,j-1,k);
            break;
        case 2:
            uidx = myWIndex(i,j,k);
            pidx0 = p_idx(i,j,k-1);
            p0 = sp.getValue(i,j,k-1);
            break;
    }

    if (!isInSystem(uidx)) {
        return 0;
    }

    if ( pidx0 == AIR && isInSystem(pidx1) ) {
        return -SYSlerp(p1, p0, uweight);
    } else if ( pidx1 == AIR && isInSystem(pidx0) ){
        return SYSlerp(p0, p1, uweight);
    }
    return 0;
}

/// =================== Implement sim_stokesSolver::classifyAndBuildIndices() ===================
template<typename T> void sim_stokesSolver<T>::classifyAndBuildIndices( const SIM_RawField * const* surf_weights,
                                                                        const SIM_RawField * const* col_weights
                                                                        ){
    initAndClassifyIndex(surf_weights, col_weights, myUIndex, FACEX);
    initAndClassifyIndex(surf_weights, col_weights, myVIndex, FACEY);
    initAndClassifyIndex(surf_weights, col_weights, myWIndex, FACEZ);

    // the central indices depend on face indices being classified.
    // We want to avoid creating pressure samples surrounded by collision faces
    initAndClassifyIndex(surf_weights, col_weights, myCentralIndex, CENTER);
    initAndClassifyIndex(surf_weights, col_weights, myTxyIndex, EDGEXY);
    initAndClassifyIndex(surf_weights, col_weights, myTxzIndex, EDGEXZ);
    initAndClassifyIndex(surf_weights, col_weights, myTyzIndex, EDGEYZ);

    exint maxindex = 0;
    buildIndex(myCentralIndex, CENTER, maxindex);
    myNumPressureVars += maxindex;
    maxindex *= 3; // account for txx and tyy indices
    buildIndex(myTxyIndex, EDGEXY, maxindex);
    buildIndex(myTxzIndex, EDGEXZ, maxindex);
    buildIndex(myTyzIndex, EDGEYZ, maxindex);
    myNumStressVars += maxindex - myNumPressureVars;

    buildVelocityIndices(surf_weights, col_weights);
}

/// =================== Implement sim_stokesSolver::initAndClassifyIndex() ===================
// Chain `classifyAndBuildIndices() -> initAndClassifyIndex()`
template<typename T> void sim_stokesSolver<T>::initAndClassifyIndex(const SIM_RawField * const* surf_weights,
                                                                    const SIM_RawField * const* col_weights,
                                                                    SIM_RawIndexField &index,
                                                                    FieldIndex fidx
                                                                    ) {
    index.match(*surf_weights[fidx]);
    index.makeConstant(INVALIDIDX);
    index.setBorder(UT_VOXELBORDER_CONSTANT, INVALIDIDX);

    classifyIndexField(surf_weights, col_weights, index, fidx);
}

/// =================== Implement sim_stokesSolver::buildIndex() ===================
// Chain `classifyAndBuildIndices() -> buildIndex()`
template<typename T> void sim_stokesSolver<T>::buildIndex(  SIM_RawIndexField &index,
                                                            FieldIndex fidx,
                                                            exint &maxindex
                                                            ){
    UT_VoxelArrayIteratorI vit(index.fieldNC());
    UT_VoxelTileIteratorI vitt;
    for (vit.rewind(); !vit.atEnd(); vit.advanceTile()) {
        if ( vit.isTileConstant() && !isInSystem(vit.getValue()) ) {
            continue;
        }

        vitt.setTile(vit);
        for (vitt.rewind(); !vitt.atEnd(); vitt.advance()) {
            if ( isInSystem(vitt.getValue()) ) {
                vitt.setValue(maxindex++);
            }
        }
    }
}


/// =================== Implement sim_stokesSolver::buildVelocityIndices() ===================
// Chain `classifyAndBuildIndices() -> buildVelocityIndices()`
// Additional indices for decoupled systems
template<typename T> void sim_stokesSolver<T>::buildVelocityIndices(const SIM_RawField * const* surf_weights,
                                                                    const SIM_RawField * const* col_weights
                                                                    ){
    // Velocity indices start from 0 as they are local to their block because they
    // are only used in the blockwise system builder
    exint maxindex = 0;
    buildIndex(myUIndex, FACEX, maxindex);
    buildIndex(myVIndex, FACEY, maxindex);
    buildIndex(myWIndex, FACEZ, maxindex);

    myCollisionIndex = maxindex;

    // build Collision Velocity indices
    buildCollisionIndex(myUIndex, FACEX, maxindex);
    buildCollisionIndex(myVIndex, FACEY, maxindex);
    buildCollisionIndex(myWIndex, FACEZ, maxindex);
    myNumVelocityVars += maxindex;
}

/// =================== Implement sim_stokesSolver::classifyIndexFieldPartial() ===================
// Chain `classifyAndBuildIndices() -> initAndClassifyIndex() -> classifyIndexField() -> classifyIndexFieldPartial()`
template<typename T> void sim_stokesSolver<T>::classifyIndexFieldPartial(   const SIM_RawField * const* surf_weights,
                                                                            const SIM_RawField * const* col_weights,
                                                                            SIM_RawIndexField &index,
                                                                            FieldIndex fidx,
                                                                            const UT_JobInfo &info
                                                                            ){
    UT_VoxelArrayIteratorI vit(index.fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);
    for (vit.rewind(); !vit.atEnd(); vit.advance()) {
        int i = vit.x(), j = vit.y(), k = vit.z();
        auto type = solveType(surf_weights, col_weights, i, j, k, fidx);
        if ( type == INVALIDIDX ) {
            continue; // already set to INVALIDIDX
        }
        vit.setValue(type);
    }
}

/// =================== Implement sim_stokesSolver::buildCollisionIndex() ===================
// Chain `classifyAndBuildIndices() -> buildVelocityIndices() -> buildCollisionIndex()`
// PRE: assume indices have already been classified
template<typename T> void sim_stokesSolver<T>::buildCollisionIndex( SIM_RawIndexField &index,
                                                                    FieldIndex fidx,
                                                                    exint &maxindex
                                                                    ){
    UT_VoxelArrayIteratorI vit(index.fieldNC());
    UT_VoxelTileIteratorI vitt;
    for (vit.rewind(); !vit.atEnd(); vit.advanceTile()) {
        if ( vit.isTileConstant() && !isCollision(vit.getValue()) ) {
            continue;
        }
        vitt.setTile(vit);
        for (vitt.rewind(); !vitt.atEnd(); vitt.advance()) {
            if ( isCollision(vitt.getValue()) ) {
                vitt.setValue(maxindex++);
            }
        }
    }
}

/// =================== Implement sim_stokesSolver::solveType() ===================
// Chain `classifyAndBuildIndices() -> initAndClassifyIndex() -> classifyIndexField() -> classifyIndexFieldPartial() -> solveType()`
template<typename T> SolveType sim_stokesSolver<T>::solveType(  const SIM_RawField * const* surf_weights,
                                                                const SIM_RawField * const* col_weights,
                                                                int i, int j, int k,
                                                                FieldIndex fidx
                                                                ) const{
    //const UT_VoxelArrayF &c_vol_liquid = *surf_weights[0]->field();
    const UT_VoxelArrayF &ez_vol_liquid = *surf_weights[1]->field();
    const UT_VoxelArrayF &ey_vol_liquid = *surf_weights[2]->field();
    const UT_VoxelArrayF &ex_vol_liquid = *surf_weights[3]->field();
    const UT_VoxelArrayF &u_vol_liquid = *surf_weights[4]->field();
    const UT_VoxelArrayF &v_vol_liquid = *surf_weights[5]->field();
    const UT_VoxelArrayF &w_vol_liquid = *surf_weights[6]->field();

    const UT_VoxelArrayF &c_vol_fluid = *col_weights[0]->field();
    const UT_VoxelArrayF &ez_vol_fluid = *col_weights[1]->field();
    const UT_VoxelArrayF &ey_vol_fluid = *col_weights[2]->field();
    const UT_VoxelArrayF &ex_vol_fluid = *col_weights[3]->field();
    const UT_VoxelArrayF &u_vol_fluid = *col_weights[4]->field();
    const UT_VoxelArrayF &v_vol_fluid = *col_weights[5]->field();
    const UT_VoxelArrayF &w_vol_fluid = *col_weights[6]->field();

    bool insystem = false;
    // check
    switch ( fidx ) {
        case FACEX:
            insystem = ( !u_oob(i,j,k) && u_vol_fluid(i,j,k) && u_vol_liquid(i,j,k) );
            break;

        case FACEY:
            insystem = ( !v_oob(i,j,k) && v_vol_fluid(i,j,k) && v_vol_liquid(i,j,k) );
            break;

        case FACEZ:
            insystem = ( !w_oob(i,j,k) && w_vol_fluid(i,j,k) && w_vol_liquid(i,j,k) );
            break;

        case CENTER:
            insystem =
                !(isCollision(myUIndex(i,j,k)) &&
                isCollision(myUIndex(i+1,j,k)) &&
                isCollision(myVIndex(i,j,k)) &&
                isCollision(myVIndex(i,j+1,k)) &&
                isCollision(myWIndex(i,j,k)) &&
                isCollision(myWIndex(i,j,k+1)) ) && // cell surrounded by walls
                ( !c_oob(i,j,k) && c_vol_fluid(i,j,k) );
            break;

        case EDGEXY:
            insystem = ( !txy_oob(i,j,k) && ez_vol_liquid(i,j,k) && ez_vol_fluid(i,j,k) );
            break;

        case EDGEXZ:
            insystem = ( !txz_oob(i,j,k) && ey_vol_liquid(i,j,k) && ey_vol_fluid(i,j,k) );
            break;

        case EDGEYZ:
            insystem = ( !tyz_oob(i,j,k) && ex_vol_liquid(i,j,k) && ex_vol_fluid(i,j,k) );
            break;

        default:
            assert(0); // unknown FIELDIDX
            break;
    }

    // not in the linear system
    if (!insystem) {
        return INVALIDIDX;
    }

    switch ( fidx ) {
        // classify collision boundary
        case FACEX:
            if ( u_oob(i+1,j,k) || u_oob(i-1,j,k) ) {
                //std::cerr << "uc at " << i << " " << j << " " << k << std::endl;
                return COLLISION;
            }
            if (  u_vol_fluid(i,j,k) < 0.5 || c_oob(i,j,k) || !c_vol_fluid(i,j,k) || c_oob(i-1,j,k) || !c_vol_fluid(i-1,j,k) ||
                !ey_vol_fluid(i,j,k) || !ey_vol_fluid(i,j,k+1) ||
                !ez_vol_fluid(i,j,k) || !ez_vol_fluid(i,j+1,k) )
                return COLLISION;
            break;

        case FACEY:
            if ( v_oob(i,j+1,k) || v_oob(i,j-1,k) ) {
                //std::cerr << "vc at " << i << " " << j << " " << k << std::endl;
                return COLLISION;
            }
            //if ( v_vol_fluid(i,j,k) < 0.5 )
            if (  v_vol_fluid(i,j,k) < 0.5 || c_oob(i,j,k) || !c_vol_fluid(i,j,k) || c_oob(i,j-1,k) || !c_vol_fluid(i,j-1,k) ||
                !ex_vol_fluid(i,j,k) || !ex_vol_fluid(i,j,k+1) ||
                !ez_vol_fluid(i,j,k) || !ez_vol_fluid(i+1,j,k) )
                return COLLISION;
            break;

        case FACEZ:
            if ( w_oob(i,j,k+1) || w_oob(i,j,k-1) ) {
                //std::cerr << "wc at " << i << " " << j << " " << k << std::endl;
                return COLLISION;
            }
            //if ( w_vol_fluid(i,j,k) < 0.5 ) return COLLISION;
            if ( w_vol_fluid(i,j,k) < 0.5 || c_oob(i,j,k) || !c_vol_fluid(i,j,k) || c_oob(i,j,k-1) || !c_vol_fluid(i,j,k-1) ||
                !ex_vol_fluid(i,j,k) || !ex_vol_fluid(i,j+1,k) ||
                !ey_vol_fluid(i,j,k) || !ey_vol_fluid(i+1,j,k) )
                return COLLISION;
            break;

        default:
            break;
    }

    switch ( fidx ) {
        case EDGEXY:
            insystem =
                u_vol_liquid(i,j,k) && !u_oob(i,j-1,k) && u_vol_liquid(i,j-1,k) &&
                v_vol_liquid(i,j,k) && !v_oob(i-1,j,k) && v_vol_liquid(i-1,j,k);
            break;
        case EDGEXZ:
            insystem =
                u_vol_liquid(i,j,k) && !u_oob(i,j,k-1) && u_vol_liquid(i,j,k-1) &&
                w_vol_liquid(i,j,k) && !w_oob(i-1,j,k) && w_vol_liquid(i-1,j,k);
            break;
        case EDGEYZ:
            insystem =
                v_vol_liquid(i,j,k) && !v_oob(i,j,k-1) && v_vol_liquid(i,j,k-1) &&
                w_vol_liquid(i,j,k) && !w_oob(i,j-1,k) && w_vol_liquid(i,j-1,k);
            break;
        case FACEX:
            assert(!c_oob(i-1,j,k) && c_vol_fluid(i-1,j,k));
            insystem =
                ( c_vol_fluid(i,j,k) &&
                ez_vol_fluid(i,j,k) && !txy_oob(i,j+1,k) && ez_vol_fluid(i,j+1,k) &&
                ey_vol_fluid(i,j,k) && !txz_oob(i,j,k+1) && ey_vol_fluid(i,j,k+1) );
            break;

        case FACEY:
            assert(!c_oob(i,j-1,k) && c_vol_fluid(i,j-1,k));
            insystem =
                ( c_vol_fluid(i,j,k) &&
                ez_vol_fluid(i,j,k) && !txy_oob(i+1,j,k) && ez_vol_fluid(i+1,j,k) &&
                ex_vol_fluid(i,j,k) && !tyz_oob(i,j,k+1) && ex_vol_fluid(i,j,k+1) );
            break;

        case FACEZ:
            assert(!c_oob(i,j,k-1) && c_vol_fluid(i,j,k-1));
            insystem =
                ( c_vol_fluid(i,j,k) &&
                ey_vol_fluid(i,j,k) && !txz_oob(i+1,j,k) && ey_vol_fluid(i+1,j,k) &&
                ex_vol_fluid(i,j,k) && !tyz_oob(i,j+1,k) && ex_vol_fluid(i,j+1,k) );
            break;

        case CENTER:
            insystem =
                (!u_oob(i+1,j,k) && u_vol_liquid(i+1,j,k)) && (!u_oob(i,j,k) && u_vol_liquid(i,j,k)) &&
                (!v_oob(i,j+1,k) && v_vol_liquid(i,j+1,k)) && (!v_oob(i,j,k) && v_vol_liquid(i,j,k)) &&
                (!w_oob(i,j,k+1) && w_vol_liquid(i,j,k+1)) && (!w_oob(i,j,k) && w_vol_liquid(i,j,k));
            break;

        default:
            break;
    }

    if ( fidx == CENTER || fidx == EDGEXY || fidx == EDGEXZ || fidx == EDGEYZ ) {
        return insystem ? SOLVED : AIR; // needed for surface tension (doesn't get an index)
    } else{
        return insystem ? SOLVED : INVALIDIDX; // needed for moving boundaries (gets an index in blockwise code)
    }
}


// minimum allowed surface weight
static const fpreal MINWEIGHT = 0.1;

/// Volume Fraction is the percentage of surface SDF area occupied by each voxel in the weight field.
static void simEstimateVolumeFractions( const SIM_RawField*  surfaceField,
                                        bool                 is_surfaceConstant,
                                        SIM_FieldSample      samplingPos,
                                        int                  nsamples,
                                        bool                 invert,
                                        SIM_RawField&        weights
                                        ) {
    UT_Vector3 size = surfaceField->getSize();
    UT_Vector3 orig = surfaceField->getOrig();
    int xres, yres, zres;
    surfaceField->getVoxelRes(xres, yres, zres);
    weights.init(samplingPos, orig, size, xres, yres, zres);

    // Compress Field if Field is Constant
    if (is_surfaceConstant){
        weights.makeConstant(1);
    }else{
        // Compute fractional volume weights representing the amount the voxel surrounding each sample point is inside the provided SDF.
        weights.computeSDFWeightsSampled(surfaceField, nsamples, invert, MINWEIGHT);
    }
}
