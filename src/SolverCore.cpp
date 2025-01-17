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

            // building the system is involved, so we need a temporary datastructure to
            // handle terms added in the same place.
            struct RowEntry {
                RowEntry(int col, T val) : col(col), val(val) { }
                ~RowEntry() { }
                bool operator<(const RowEntry& other) const { return col < other.col; } // column comparator
                int col;
                T val;
            };

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

            void addUTerm(  int row_index, int i, int j, int k, float sign, float outer_liquid, float outer_fluid,
                            const sim_buildSystemParms& parms,
                            UT_VoxelProbeAverage<float,-1,0,0>& rhox,
                            UT_Array<RowEntry>& rowentries,
                            VectorType& b
                            ) const;
            void addVTerm(  int row_index, int i, int j, int k, float sign, float outer_liquid, float outer_fluid,
                            const sim_buildSystemParms& parms,
                            UT_VoxelProbeAverage<float,0,-1,0>& rhoy,
                            UT_Array<RowEntry>& rowentries,
                            VectorType& b
                            ) const;
            void addWTerm(  int row_index, int i, int j, int k, float sign, float outer_liquid, float outer_fluid,
                            const sim_buildSystemParms& parms,
                            UT_VoxelProbeAverage<float,0,0,-1>& rhoz,
                            UT_Array<RowEntry>& rowentries,
                            VectorType& b
                            ) const;

            // we need this function to combine entries with the same column index
            void appendRowEntries(MatrixType& A, int row_index, UT_Array<RowEntry>& rowentries) const;

            // pressure block:
            exint p_blk_idx(int i, int j, int k) const { return p_idx(i,j,k); }

            // stress tensor block:
            exint txx_blk_idx(int i, int j, int k) const { return txx_idx(i,j,k) - myNumPressureVars; }
            exint tyy_blk_idx(int i, int j, int k) const { return tyy_idx(i,j,k) - myNumPressureVars; }
            exint tyz_blk_idx(int i, int j, int k) const { return tyz_idx(i,j,k) - myNumPressureVars; }
            exint txz_blk_idx(int i, int j, int k) const { return txz_idx(i,j,k) - myNumPressureVars; }
            exint txy_blk_idx(int i, int j, int k) const { return txy_idx(i,j,k) - myNumPressureVars; }

            exint tzz_blk_idx(int i, int j, int k) const {
                return c_oob(i,j,k) ? exint(INVALIDIDX) : myCentralIndex(i,j,k) + myNumStressVars - myNumPressureVars;
            }

            // velocity block: ( this is not in the final system, but used to build
            // intermediate operators, like deformation rate operator, and gradient
            // operator )
            exint u_blk_idx(int i, int j, int k) const {
                return u_oob(i,j,k) ? exint(INVALIDIDX) : myUIndex(i,j,k);
            }
            exint v_blk_idx(int i, int j, int k) const {
                return v_oob(i,j,k) ? exint(INVALIDIDX) : myVIndex(i,j,k);
            }
            exint w_blk_idx(int i, int j, int k) const {
                return w_oob(i,j,k) ? exint(INVALIDIDX) : myWIndex(i,j,k);
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
    std::cerr << " Solve! " << std::endl;
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
    std::cerr << " Solve Stokes! " << std::endl;
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
    std::cerr << " Build System! " << std::endl;
    std::cerr << " Initialize b by the Zero " << std::endl;
    b.zero();

    std::cerr << " Add Center Terms " << std::endl;
    addCenterTerms(A,b,parms);

    std::cerr << " Add Txy Terms " << std::endl;
    addTxyTerms(A,b,parms);

    std::cerr << " Add Txz Terms " << std::endl;
    addTxzTerms(A,b,parms);

    std::cerr << " Add Tyz Terms " << std::endl;
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
    std::cerr << " Solve System! " << std::endl;
    b.testForNan();
    auto system_size = b.length();
    if(!(b.length() == A.getNumRows())) {
        return FAILED;
    }
    if ( !system_size ) {
        return NOCHANGE;
    }

    UT_PerfMonAutoSolveEvent event(&mySolver, "Solve Stokes Linear System");
    std::cerr << " Solve Stokes Linear System! " << std::endl;
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
    std::cerr << " Update Velocities! " << std::endl;
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
        std::cerr << " Sampling Axis is 1! " << std::endl;
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
        std::cerr << " Sampling Axis is 2! " << std::endl;
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
        std::cerr << " Sampling Axis is Unknown Dimension(3)! " << std::endl;
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
                                                                        ) {
    std::cerr << " Classify and Build Indices! " << std::endl;
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
    std::cerr << " Init and Classify Index! " << std::endl;
    std::cerr << "surf_weights[" << fidx << "] = " << surf_weights[fidx] << std::endl;
    std::cerr << "fidx = " << fidx << std::endl;
    index.match(*surf_weights[fidx]);
    index.makeConstant(INVALIDIDX);
    index.setBorder(UT_VOXELBORDER_CONSTANT, INVALIDIDX);

    std::cerr << " Call classifyIndexField! " << std::endl;
    classifyIndexField(surf_weights, col_weights, index, fidx);
}

/// =================== Implement sim_stokesSolver::buildIndex() ===================
// Chain `classifyAndBuildIndices() -> buildIndex()`
template<typename T> void sim_stokesSolver<T>::buildIndex(  SIM_RawIndexField &index,
                                                            FieldIndex fidx,
                                                            exint &maxindex
                                                            ){
    std::cerr << " Build Index! " << std::endl;
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
    std::cerr << " Build Velocity Indices! " << std::endl;
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
    std::cerr << " Classify Index Field Partial! " << std::endl;
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
    std::cerr << " Build Collision Index! " << std::endl;
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
    //std::cerr << " Solve Type! " << std::endl;
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


/// =================== Implement sim_stokesSolver::addCenterTermsPartial() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addCenterTerms() -> addCenterTermsPartial()`
template<typename T> void sim_stokesSolver<T>::addCenterTermsPartial(   MatrixType &A,
                                                                        VectorType &b,
                                                                        const sim_buildSystemParms& parms,
                                                                        const UT_JobInfo& info
                                                                        ) const {
    // Setup density probes
    UT_VoxelProbeAverage<float,-1,0,0> rho_x;
    UT_VoxelProbeAverage<float,0,-1,0> rho_y;
    UT_VoxelProbeAverage<float,0,0,-1> rho_z;
    rho_x.setArray(&parms.density);
    rho_y.setArray(&parms.density);
    rho_z.setArray(&parms.density);

    auto min_visc = mySolver.getMinViscosity();

    //rhs << Bp*uold - dx*WLp*(G.transpose()*WFu - WFp*G.transpose())*ubc + WLp*G.transpose()*WFu*ust,
    //       Bt*uold - dx*WLt*(D*WFu - WFt*D)*ubc + WLt*D*WFu*ust;
    //
    UT_Array<RowEntry> rowentries;
    UT_VoxelArrayIteratorI vit;
    UT_VoxelTileIteratorI vitt;
    vit.setConstArray(myCentralIndex.field());
    vit.splitByTile(info);
    for ( vit.rewind(); !vit.atEnd(); vit.advanceTile() ) {
        if ( vit.isTileConstant() && !isInSystem(vit.getValue()) ) {
            continue;
        }

        vitt.setTile(vit);

        for ( vitt.rewind(); !vitt.atEnd(); vitt.advance() ) {
            if (!isInSystem(vitt.getValue())) {
                continue;
            }
            int i = vitt.x(), j = vitt.y(), k = vitt.z();

            auto cfw = parms.c_vol_fluid(i,j,k); // central fluid volume weight
            auto clw = parms.c_vol_liquid(i,j,k); // central liquid volume weight

            // du/dx + dv/dy + dw/dz = 0
            int row_index = p_idx(i,j,k);
            rowentries.clear(); // reset entries per row

            b(row_index) = 0;

            addUTerm(row_index, i+1, j, k, +1, clw, cfw, parms, rho_x, rowentries, b);
            addUTerm(row_index, i,   j, k, -1, clw, cfw, parms, rho_x, rowentries, b);

            addVTerm(row_index, i, j+1, k, +1, clw, cfw, parms, rho_y, rowentries, b);
            addVTerm(row_index, i, j,   k, -1, clw, cfw, parms, rho_y, rowentries, b);

            addWTerm(row_index, i, j, k+1, +1, clw, cfw, parms, rho_z, rowentries, b);
            addWTerm(row_index, i, j, k,   -1, clw, cfw, parms, rho_z, rowentries, b);

            appendRowEntries(A, row_index, rowentries);

            // txx + 0.5*tyy - du/dx + dw/dz = 0
            row_index = txx_idx(i,j,k);
            rowentries.clear();

            auto visc = parms.viscosity(i,j,k);
            auto factor = visc < min_visc ? 0 : dx*dx/visc;
            auto diag = clw * cfw * factor;

            rowentries.emplace_back(row_index, diag);
            rowentries.emplace_back(tyy_idx(i,j,k), 0.5f*diag);
            b(row_index) = 0;

            addUTerm(row_index, i+1,j, k,  -1, clw, cfw, parms, rho_x, rowentries, b);
            addUTerm(row_index, i,  j, k,  +1, clw, cfw, parms, rho_x, rowentries, b);

            addWTerm(row_index, i, j, k+1, +1, clw, cfw, parms, rho_z, rowentries, b);
            addWTerm(row_index, i, j, k,   -1, clw, cfw, parms, rho_z, rowentries, b);

            appendRowEntries(A, row_index, rowentries);

            // tyy + 0.5*txx - dv/dy + dw/dz = 0
            row_index = tyy_idx(i,j,k);
            rowentries.clear();

            rowentries.emplace_back(row_index, diag);
            rowentries.emplace_back(txx_idx(i,j,k), 0.5*diag);
            b(row_index) = 0;

            addVTerm(row_index, i, j+1, k, -1, clw, cfw, parms, rho_y, rowentries, b);
            addVTerm(row_index, i, j,   k, +1, clw, cfw, parms, rho_y, rowentries, b);

            addWTerm(row_index, i, j, k+1, +1, clw, cfw, parms, rho_z, rowentries, b);
            addWTerm(row_index, i, j, k,   -1, clw, cfw, parms, rho_z, rowentries, b);

            appendRowEntries(A, row_index, rowentries);
        }
    }
}


/// =================== Implement sim_stokesSolver::addTxyTermsPartial() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addTxyTerms() -> addTxyTermsPartial()`
template<typename T> void sim_stokesSolver<T>::addTxyTermsPartial(  MatrixType &A,
                                                                    VectorType &b,
                                                                    const sim_buildSystemParms& parms,
                                                                    const UT_JobInfo& info) const {
    UT_VoxelProbeAverage<float, -1, -1, 0> visc_xy;
    visc_xy.setArray(&parms.viscosity);

    auto min_visc = mySolver.getMinViscosity();

    UT_VoxelProbeAverage<float,-1,0,0> rho_x;
    UT_VoxelProbeAverage<float,0,-1,0> rho_y;
    rho_x.setArray(&parms.density);
    rho_y.setArray(&parms.density);

    UT_Array<RowEntry> rowentries;
    //txy - du/dy - dv/dx = 0
    UT_VoxelArrayIteratorI vit;
    UT_VoxelTileIteratorI vitt;
    vit.setConstArray(myTxyIndex.field());
    vit.splitByTile(info);
    for ( vit.rewind(); !vit.atEnd(); vit.advanceTile() ) {
        if ( vit.isTileConstant() && !isInSystem(vit.getValue()) ) {
            continue;
        }

        vitt.setTile(vit);

        for ( vitt.rewind(); !vitt.atEnd(); vitt.advance() ) {
            if (!isInSystem(vitt.getValue())) {
                continue;
            }

            int i = vitt.x(), j = vitt.y(), k = vitt.z();
            int row_index = txy_idx(i,j,k);
            rowentries.clear();

            visc_xy.setIndex(vitt);
            auto visc = visc_xy.getValue();
            auto factor = visc < min_visc ? 0.0 : dx*dx/visc;
            auto lw = parms.ez_vol_liquid(i,j,k);
            auto fw = parms.ez_vol_fluid(i,j,k);
            rowentries.emplace_back(row_index, lw * fw * factor);
            b(row_index) = 0;

            addUTerm(row_index, i, j,   k, -1, lw, fw, parms, rho_x, rowentries, b);
            addUTerm(row_index, i, j-1, k, +1, lw, fw, parms, rho_x, rowentries, b);

            addVTerm(row_index, i,   j, k, -1, lw, fw, parms, rho_y, rowentries, b);
            addVTerm(row_index, i-1, j, k, +1, lw, fw, parms, rho_y, rowentries, b);

            appendRowEntries(A, row_index, rowentries);
        }
    }
}


/// =================== Implement sim_stokesSolver::addTxzTermsPartial() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addTxzTerms() -> addTxzTermsPartial()`
template<typename T> void sim_stokesSolver<T>::addTxzTermsPartial(  MatrixType &A,
                                                                    VectorType &b,
                                                                    const sim_buildSystemParms& parms,
                                                                    const UT_JobInfo& info
                                                                    ) const {
    UT_VoxelProbeAverage<float, -1, 0, -1> visc_xz;
    visc_xz.setArray(&parms.viscosity);

    auto min_visc = mySolver.getMinViscosity();

    UT_VoxelProbeAverage<float,-1,0,0> rho_x;
    UT_VoxelProbeAverage<float,0,0,-1> rho_z;
    rho_x.setArray(&parms.density);
    rho_z.setArray(&parms.density);

    UT_Array<RowEntry> rowentries;
    //txz - du/dz - dw/dx = 0
    UT_VoxelArrayIteratorI vit;
    UT_VoxelTileIteratorI vitt;
    vit.setConstArray(myTxzIndex.field());
    vit.splitByTile(info);
    for ( vit.rewind(); !vit.atEnd(); vit.advanceTile() ) {
        if ( vit.isTileConstant() && !isInSystem(vit.getValue()) ) {
            continue;
        }


        vitt.setTile(vit);
        for ( vitt.rewind(); !vitt.atEnd(); vitt.advance() ) {
            if (!isInSystem(vitt.getValue())) {
                continue;
            }

            int i = vitt.x(), j = vitt.y(), k = vitt.z();
            int row_index = txz_idx(i,j,k);
            rowentries.clear();

            visc_xz.setIndex(vitt);
            auto visc = visc_xz.getValue();
            auto factor = visc < min_visc ? 0.0 : dx*dx/visc;
            auto lw = parms.ey_vol_liquid(i,j,k);
            auto fw = parms.ey_vol_fluid(i,j,k);
            rowentries.emplace_back(row_index, lw * fw * factor);
            b(row_index) = 0;

            addUTerm(row_index, i, j, k,   -1, lw, fw, parms, rho_x, rowentries, b);
            addUTerm(row_index, i, j, k-1, +1, lw, fw, parms, rho_x, rowentries, b);

            addWTerm(row_index, i,   j, k, -1, lw, fw, parms, rho_z, rowentries, b);
            addWTerm(row_index, i-1, j, k, +1, lw, fw, parms, rho_z, rowentries, b);

            appendRowEntries(A, row_index, rowentries);
        }
    }
}


/// =================== Implement sim_stokesSolver::addTyzTermsPartial() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addTyzTerms() -> addTyzTermsPartial()`
template<typename T> void sim_stokesSolver<T>::addTyzTermsPartial(  MatrixType &A,
                                                                    VectorType &b,
                                                                    const sim_buildSystemParms& parms,
                                                                    const UT_JobInfo& info
                                                                    ) const {
    UT_VoxelProbeAverage<float, 0, -1, -1> visc_yz;
    visc_yz.setArray(&parms.viscosity);

    auto min_visc = mySolver.getMinViscosity();

    UT_VoxelProbeAverage<float,0,-1,0> rho_y;
    UT_VoxelProbeAverage<float,0,0,-1> rho_z;
    rho_y.setArray(&parms.density);
    rho_z.setArray(&parms.density);

    UT_Array<RowEntry> rowentries;
    //tyz = dv/dz + dw/dy
    UT_VoxelArrayIteratorI vit;
    UT_VoxelTileIteratorI vitt;
    vit.setConstArray(myTyzIndex.field());
    vit.splitByTile(info);
    for ( vit.rewind(); !vit.atEnd(); vit.advanceTile() ) {
        if ( vit.isTileConstant() && !isInSystem(vit.getValue()) ) {
            continue;
        }

        vitt.setTile(vit);
        for ( vitt.rewind(); !vitt.atEnd(); vitt.advance() ) {
            if (!isInSystem(vitt.getValue())) {
                continue;
            }

            int i = vitt.x(), j = vitt.y(), k = vitt.z();
            int row_index = tyz_idx(i,j,k);
            rowentries.clear();

            visc_yz.setIndex(vitt);
            auto visc = visc_yz.getValue();
            auto factor = visc < min_visc ? 0.0 : dx*dx/visc;
            auto lw = parms.ex_vol_liquid(i,j,k);
            auto fw = parms.ex_vol_fluid(i,j,k);
            rowentries.emplace_back(row_index, lw * fw * factor);
            b(row_index) = 0;

            addVTerm(row_index, i, j,   k,   -1, lw, fw, parms, rho_y, rowentries, b);
            addVTerm(row_index, i, j,   k-1, +1, lw, fw, parms, rho_y, rowentries, b);

            addWTerm(row_index, i, j,   k,   -1, lw, fw, parms, rho_z, rowentries, b);
            addWTerm(row_index, i, j-1, k,   +1, lw, fw, parms, rho_z, rowentries, b);

            appendRowEntries(A, row_index, rowentries);
        }
    }
}


/// =================== Implement sim_stokesSolver::addUTerm() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addCenterTerms() -> addCenterTermsPartial() -> addUTerm()`
template<typename T> void sim_stokesSolver<T>::addUTerm(int row_index, int i, int j, int k, float sign, float outer_liquid, float outer_fluid,
                                                        const sim_buildSystemParms& parms,
                                                        UT_VoxelProbeAverage<float,-1,0,0>& rhox,
                                                        UT_Array<RowEntry>& rowentries,
                                                        VectorType& b
                                                        ) const {
    auto solid_vel = parms.u_solid.getValue(i,j,k);
    auto vel_fw = parms.u_vol_fluid(i,j,k); // x-face fluid volume weight
    auto vel_lw = parms.u_vol_liquid(i,j,k); // x-face liquid volume weight

    int idx = myUIndex(i,j,k);
    if (isCollision(idx)){
        b(row_index) -= sign * outer_liquid * vel_fw * solid_vel * dx;
    } else if (isInSystem(idx)) {
        rhox.setIndex(i,j,k);
        auto rho = SYSclamp(rhox.getValue(), parms.minrho, parms.maxrho);
        double factor = sign * dt * outer_liquid * vel_fw / (rho * vel_lw);

        // clang-format off
        //-dp/dx
        if(isInSystem(p_idx(i,   j,   k)))     rowentries.emplace_back(p_idx(i,  j,  k),   -factor * parms.c_vol_liquid(i,  j,  k));
        if(isInSystem(p_idx(i-1, j,   k)))     rowentries.emplace_back(p_idx(i-1,j,  k),   +factor * parms.c_vol_liquid(i-1,j,  k));

        //dtxx/dx
        if(isInSystem(txx_idx(i,  j,  k)))   rowentries.emplace_back(txx_idx(i,  j,  k),   +factor * parms.c_vol_liquid(i,  j,  k));
        if(isInSystem(txx_idx(i-1,j,  k)))   rowentries.emplace_back(txx_idx(i-1,j,  k),   -factor * parms.c_vol_liquid(i-1,j,  k));

        //dtxy/dy
        if(isInSystem(txy_idx(i,  j+1,k)))   rowentries.emplace_back(txy_idx(i,  j+1,k),   +factor * parms.ez_vol_liquid(i, j+1,k));
        if(isInSystem(txy_idx(i,  j,  k)))   rowentries.emplace_back(txy_idx(i,  j,  k),   -factor * parms.ez_vol_liquid(i, j,  k));

        //dtxz/dz
        if(isInSystem(txz_idx(i,  j,  k+1))) rowentries.emplace_back(txz_idx(i,  j,  k+1), +factor * parms.ey_vol_liquid(i, j,  k+1));
        if(isInSystem(txz_idx(i,  j,  k)))   rowentries.emplace_back(txz_idx(i,  j,  k),   -factor * parms.ey_vol_liquid(i, j,  k));
        // clang-format on

        //u*
        b(row_index) -= sign * outer_liquid * vel_fw * parms.u(i,j,k) * dx;

        auto gfp = ghostFluidSurfaceTensionPressure<0>(i,j,k, vel_lw, parms.surfpres);
        b(row_index) += factor * gfp;
    }

    b(row_index) += sign * outer_liquid * vel_fw      * solid_vel * dx;
    b(row_index) -= sign * outer_liquid * outer_fluid * solid_vel * dx;
}


/// =================== Implement sim_stokesSolver::addVTerm() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addCenterTerms() -> addCenterTermsPartial() -> addVTerm()`
template<typename T> void sim_stokesSolver<T>::addVTerm(int row_index, int i, int j, int k, float sign, float outer_liquid, float outer_fluid,
                                                        const sim_buildSystemParms& parms,
                                                        UT_VoxelProbeAverage<float,0,-1,0>& rhoy,
                                                        UT_Array<RowEntry>& rowentries,
                                                        VectorType& b
                                                        ) const {
    auto solid_vel = parms.v_solid.getValue(i,j,k);
    auto vel_fw = parms.v_vol_fluid(i,j,k); // y-face fluid volume weight
    auto vel_lw = parms.v_vol_liquid(i,j,k); // y-face liquid volume weight

    int idx = myVIndex(i,j,k);
    if (isCollision(idx)) {
        b(row_index) -=  sign * outer_liquid * vel_fw * solid_vel * dx;
    } else if (isInSystem(idx)) {
        rhoy.setIndex(i,j,k);
        auto rho = SYSclamp(rhoy.getValue(), parms.minrho, parms.maxrho);
        double factor = sign * dt * outer_liquid * vel_fw / (rho * vel_lw);

        // clang-format off
        //-dp/dy
        if(isInSystem(p_idx(i,   j,   k)))   rowentries.emplace_back(  p_idx(i,  j,  k),   -factor * parms.c_vol_liquid(i,   j,  k));
        if(isInSystem(p_idx(i,   j-1, k)))   rowentries.emplace_back(  p_idx(i,  j-1,k),   +factor * parms.c_vol_liquid(i,   j-1,k));

        //dtxy/dx
        if(isInSystem(txy_idx(i+1,j,  k)))   rowentries.emplace_back(txy_idx(i+1,j,  k),   +factor * parms.ez_vol_liquid(i+1,j,  k));
        if(isInSystem(txy_idx(i,  j,  k)))   rowentries.emplace_back(txy_idx(i,  j,  k),   -factor * parms.ez_vol_liquid(i,  j,  k));

        //dtyy/dy
        if(isInSystem(tyy_idx(i,  j,  k)))   rowentries.emplace_back(tyy_idx(i,  j,  k),   +factor * parms.c_vol_liquid(i,   j,  k));
        if(isInSystem(tyy_idx(i,  j-1,k)))   rowentries.emplace_back(tyy_idx(i,  j-1,k),   -factor * parms.c_vol_liquid(i,   j-1,k));

        //dtyz/dz
        if(isInSystem(tyz_idx(i,  j,  k+1))) rowentries.emplace_back(tyz_idx(i,  j,  k+1), +factor * parms.ex_vol_liquid(i,   j, k+1));
        if(isInSystem(tyz_idx(i,  j,  k)))   rowentries.emplace_back(tyz_idx(i,  j,  k),   -factor * parms.ex_vol_liquid(i,   j, k));
        // clang-format on

        b(row_index) -= sign * outer_liquid * vel_fw * parms.v(i,j,k) * dx;

        auto gfp = ghostFluidSurfaceTensionPressure<1>(i,j,k, vel_lw, parms.surfpres);
        b(row_index) += factor * gfp;
    }

    b(row_index) += sign * outer_liquid * vel_fw      * solid_vel * dx;
    b(row_index) -= sign * outer_liquid * outer_fluid * solid_vel * dx;
}


/// =================== Implement sim_stokesSolver::addWTerm() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addCenterTerms() -> addCenterTermsPartial() -> addWTerm()`
template<typename T> void sim_stokesSolver<T>::addWTerm(int row_index, int i, int j, int k, float sign, float outer_liquid, float outer_fluid,
                                                        const sim_buildSystemParms& parms,
                                                        UT_VoxelProbeAverage<float,0,0,-1>& rhoz,
                                                        UT_Array<RowEntry>& rowentries,
                                                        VectorType& b
                                                        ) const {
    auto solid_vel = parms.w_solid.getValue(i,j,k);
    auto vel_fw = parms.w_vol_fluid(i,j,k); // y-face fluid volume weight
    auto vel_lw = parms.w_vol_liquid(i,j,k); // y-face liquid volume weight

    int idx = myWIndex(i,j,k);
    if (isCollision(idx)) {
        b(row_index) -= sign * outer_liquid * vel_fw * solid_vel * dx;
    } else if (isInSystem(idx)) {
        rhoz.setIndex(i,j,k);
        auto rho = SYSclamp(rhoz.getValue(), parms.minrho, parms.maxrho);
        double factor = sign * dt * outer_liquid * vel_fw / (rho * vel_lw);

        // clang-format off
        //-dpdz
        if(isInSystem(p_idx(i,   j,   k)))   rowentries.emplace_back(  p_idx(i,  j,  k),   -factor * parms.c_vol_liquid(i,   j,  k));
        if(isInSystem(p_idx(i,   j,   k-1))) rowentries.emplace_back(  p_idx(i,  j,  k-1), +factor * parms.c_vol_liquid(i,   j,  k-1));

        //dtxz/dx
        if(isInSystem(txz_idx(i+1,j,  k)))   rowentries.emplace_back(txz_idx(i+1,j,  k),   +factor * parms.ey_vol_liquid(i+1,j,  k));
        if(isInSystem(txz_idx(i,  j,  k)))   rowentries.emplace_back(txz_idx(i,  j,  k),   -factor * parms.ey_vol_liquid(i,  j,  k));

        //dtyz/dy
        if(isInSystem(tyz_idx(i,  j+1,k)))   rowentries.emplace_back(tyz_idx(i,  j+1,k),   +factor * parms.ex_vol_liquid(i,  j+1,k));
        if(isInSystem(tyz_idx(i,  j,  k)))   rowentries.emplace_back(tyz_idx(i,  j,  k),   -factor * parms.ex_vol_liquid(i,  j,  k));

        //dtzz/dz -> -dtxx/dz - dtyy/dz
        if(isInSystem(txx_idx(i,  j,  k)))   rowentries.emplace_back(txx_idx(i,  j,  k),   -factor * parms.c_vol_liquid(i,   j,  k));
        if(isInSystem(txx_idx(i,  j,  k-1))) rowentries.emplace_back(txx_idx(i,  j,  k-1), +factor * parms.c_vol_liquid(i,   j,  k-1));

        if(isInSystem(tyy_idx(i,  j,  k)))   rowentries.emplace_back(tyy_idx(i,  j,  k),   -factor * parms.c_vol_liquid(i,   j,  k));
        if(isInSystem(tyy_idx(i,  j,  k-1))) rowentries.emplace_back(tyy_idx(i,  j,  k-1), +factor * parms.c_vol_liquid(i,   j,  k-1));
        // clang-format on

        b(row_index) -= sign * outer_liquid * vel_fw * parms.w(i,j,k) * dx;

        auto gfp = ghostFluidSurfaceTensionPressure<2>(i,j,k, vel_lw, parms.surfpres);
        b(row_index) += factor * gfp;
    }

    b(row_index) += sign * outer_liquid * vel_fw      * solid_vel * dx;
    b(row_index) -= sign * outer_liquid * outer_fluid * solid_vel * dx;
}


/// =================== Implement sim_stokesSolver::appendRowEntries() ===================
// Chain `solve() -> solveStokes() -> buildSystem() -> addTxzTerms() -> addTxzTermsPartial() -> appendRowEntries()`
// we need this function to add entries with the same column index
template<typename T> void sim_stokesSolver<T>::appendRowEntries(MatrixType& A, int row_index, UT_Array<RowEntry>& rowentries) const {
    rowentries.sort(std::less<RowEntry>());
    auto it = rowentries.begin();
    if ( it == rowentries.end() ) return; // empty
    auto prev = it;
    int nz = 0;
    for ( ++it; it != rowentries.end(); ++it ) {
        if ( it->col != prev->col ) {
            A.appendRowElement(row_index, prev->col, prev->val, nz);
            prev = it;
        } else {
            prev->val += it->val;
        }
    }
    A.appendRowElement(row_index, prev->col, prev->val, nz); // handle the last element
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
