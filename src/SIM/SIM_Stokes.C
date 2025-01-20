
//#define BLOCKWISE_STOKES
//#define USE_EIGEN_SOLVER_FOR_BLOCKWISE_STOKES
//#define PRINT_ROTATING_BALL_ANGULAR_MOMENTUM

using std::tuple;
using std::make_tuple;
using std::tie;
using std::move;

namespace Stokes{

// import all enum values for convenience
using namespace Stokes;

template<typename T> class sim_stokesSolver {

    public:



        // return true if velocity index represents a collision velocity in the system
        bool isCollision(exint idx) const {
            return idx == COLLISION || idx >= myCollisionIndex;
        }







        void buildSystemBlockwise(  BlockMatrixType &matrix, BlockVectorType &rhs, BlockMatrixType &H, BlockVectorType &ust,
                                    const BlockVectorType &ustar,
                                    const SIM_RawField & surf,
                                    const SIM_RawField * const* surf_weights,
                                    const SIM_RawField * const* col_weights,
                                    const SIM_RawField &viscosity,
                                    const SIM_RawField &density,
                                    const SIM_RawField * const* solid_vel,
                                    const SIM_RawField & surf_pres
                                    ) const;
        void buildDecoupledSystem(  BlockMatrixType &,  BlockMatrixType &, BlockMatrixType &,
                                    BlockMatrixType &,  BlockMatrixType &, BlockMatrixType &,
                                    const SIM_RawField & surf,
                                    const SIM_RawField * const* surf_weights,
                                    const SIM_RawField * const* col_weights,
                                    const SIM_RawField &viscosity,
                                    const SIM_RawField &density,
                                    const SIM_RawField * const* solid_vel,
                                    const SIM_RawField & surf_pres
                                    ) const;
        void buildPressureOnlySystem(   BlockMatrixType &,  BlockMatrixType &, BlockMatrixType &,
                                        const SIM_RawField & surf,
                                        const SIM_RawField * const* surf_weights,
                                        const SIM_RawField * const* col_weights,
                                        const SIM_RawField &density,
                                        const SIM_RawField * const* solid_vel,
                                        const SIM_RawField & surf_pres
                                        ) const;
        void buildViscositySystem(  BlockMatrixType &,  BlockMatrixType &,
                                    const SIM_RawField * const* surf_weights,
                                    const SIM_RawField * const* col_weights,
                                    const SIM_RawField &viscosity,
                                    const SIM_RawField &density,
                                    const SIM_RawField * const* solid_vel
                                    ) const;
        void assembleBlockSystem(   const BlockMatrixType& WLp,
                                    const BlockMatrixType& WLuinv,
                                    const BlockMatrixType& WFu,
                                    const BlockMatrixType& WLt,
                                    const BlockMatrixType& WFt,
                                    const BlockMatrixType& G,
                                    const BlockMatrixType& D,
                                    const BlockMatrixType& Pinv,
                                    const BlockMatrixType& Minv,
                                    BlockMatrixType& Ap,
                                    BlockMatrixType& Bp,
                                    BlockMatrixType& Hp,
                                    BlockMatrixType& At,
                                    BlockMatrixType& Bt,
                                    BlockMatrixType& Ht
                                    ) const;
        void assembleStressVelocitySystem(  const BlockMatrixType& WLt,
                                            const BlockMatrixType& WLu,
                                            const BlockMatrixType& WFtinv,
                                            const BlockMatrixType& WFu,
                                            const BlockMatrixType& D,
                                            const BlockMatrixType& Pinv,
                                            const BlockMatrixType& M,
                                            BlockMatrixType &A,
                                            BlockMatrixType &B
                                            ) const;
        void removeNullSpace(const MatrixType &matrix, const VectorType &rhs) const;

    // remove zero rows and columns from A
        void pruneSystem(   const BlockMatrixType &A,
                            const BlockVectorType &b,
                            MatrixType& newA,
                            VectorType &newb,
                            UT_ExintArray& to_original
                            ) const;
        void copySystem(const BlockMatrixType &A,
                        const BlockVectorType &b,
                        MatrixType& newA,
                        VectorType &newb
                        ) const;
        SolverResult solveBlockwiseStokes(  const SIM_RawField & surf,
                                            const SIM_RawField * const* sweights,
                                            const SIM_RawField * const* cweights,
                                            const SIM_RawField & viscosity,
                                            const SIM_RawField & density,
                                            const SIM_RawField * const* solid_vel,
                                            const SIM_RawField & surf_pres,
                                            SIM_VectorField * valid,
                                            SIM_VectorField & vel
                                            ) const;


        SolverResult solveSystemEigen(  const BlockMatrixType &A,
                                        const BlockVectorType &b,
                                        BlockVectorType &x
                                        ) const;


        void updateVelocitiesBlockwise( const VecX<T> &x,
                                        const SIM_RawField * const* solid_vel,
                                        SIM_VectorField *valid,
                                        SIM_VectorField &vel
                                        ) const;

    // interpolate ghost fluid pressure at the liquid surface inside the given
    // velocity voxel


        auto buildVelocityVector(const SIM_VectorField &vel, const SIM_RawField * const* colvel) const -> BlockVectorType;
        auto buildSolidVelocityVector(const SIM_RawField * const* vel) const -> BlockVectorType;
        auto buildSurfaceTensionPressureVector(const SIM_RawField & surfp) const -> BlockVectorType;
        auto buildGhostFluidSurfaceTensionPressureVector(const SIM_RawField * const* surf_weights, const SIM_RawField & surfp) const -> BlockVectorType;
        auto buildSurfaceTensionRHSAlt(const SIM_RawField * const* surf_weights, const SIM_RawField & density, const BlockVectorType &pbc) const -> BlockVectorType;
        auto buildSurfaceTensionRHS(const SIM_RawField * const* surf_weights, const SIM_RawField & density, const BlockVectorType &ust) const -> BlockVectorType;


    private: // routine members
        // System builder helpers










        bool remove_expansion_rate_tensor() const {
            return myScheme == DECOUPLED_NOEXPANSION || myScheme == DECOUPLED_NOEXPANSION_FANCY;
        }







        // Additional index accessors provided for decoupled systems wrt their
        // corresponding block (pressure, stress and velocity blocks have independent indices)



    public:
        void buildDeformationRateOperator(BlockMatrixType& D) const;
        void buildGradientOperator(BlockMatrixType& G) const;
        void buildGhostFluidMatrix( const UT_VoxelArrayF & u_weights,
                                    const UT_VoxelArrayF & v_weights,
                                    const UT_VoxelArrayF & w_weights,
                                    BlockMatrixType& GF
                                    ) const;
        void buildSumNeighboursOperator(BlockMatrixType& N) const;
        template<bool INVERSE> void buildViscosityMatrix(const SIM_RawField & viscosity, BlockMatrixType& M) const;
        void buildDensityMatrix(const SIM_RawField & density,
                                BlockMatrixType &P
                                ) const;
        void buildPressureWeightMatrix(const UT_VoxelArrayF &c_weights, BlockMatrixType& Wp) const;
        template<bool INVERSE> void buildVelocityWeightMatrix(  const UT_VoxelArrayF &u_weights,
                                                                const UT_VoxelArrayF &v_weights,
                                                                const UT_VoxelArrayF &w_weights,
                                                                BlockMatrixType& Wu
                                                                ) const;
        template<bool INVERSE> void buildStressWeightMatrix(const UT_VoxelArrayF &c_weights,
                                                            const UT_VoxelArrayF &ex_weights,
                                                            const UT_VoxelArrayF &ey_weights,
                                                            const UT_VoxelArrayF &ez_weights,
                                                            BlockMatrixType& Wt
                                                            ) const;


        int getNumPressureVars() const { return myNumPressureVars; }
        int getNumStressVars() const { return myNumStressVars; }
        int getNumVelocityVars() const {
            return myNumVelocityVars;
        }


    private:
    // workspace triplets for use in non multithreaded functions
    mutable std::vector<Triplet<T>> triplets;
    };
} // namespace

/// Standard constructor, note that BaseClass was crated by the
/// DECLARE_DATAFACTORY and provides an easy way to chain through
/// the class hierarchy.




















template<typename T>
void
sim_stokesSolver<T>::buildDeformationRateOperator(BlockMatrixType &D) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();
  auto set_val = [&](int row, int col, T val)
  {
    if ( isInSystem(col) )
      triplets.emplace_back( row, col, val );
  };

  // cell center values: txx and tyy
  UT_VoxelArrayIteratorI vit;
  vit.setConstArray(myCentralIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if ( !isInSystem(vit.getValue()) )
      continue;

    if ( reduced_stress_tensor() )
    {
      // txx
      set_val( txx_blk_idx(i,j,k), u_blk_idx(i,j,k),   -1);
      set_val( txx_blk_idx(i,j,k), u_blk_idx(i+1,j,k),  1);
      set_val( txx_blk_idx(i,j,k), w_blk_idx(i,j,k),    1);
      set_val( txx_blk_idx(i,j,k), w_blk_idx(i,j,k+1), -1);

      // tyy
      set_val( tyy_blk_idx(i,j,k), v_blk_idx(i,j,k),   -1);
      set_val( tyy_blk_idx(i,j,k), v_blk_idx(i,j+1,k),  1);
      set_val( tyy_blk_idx(i,j,k), w_blk_idx(i,j,k),    1);
      set_val( tyy_blk_idx(i,j,k), w_blk_idx(i,j,k+1), -1);
    }
    else
    {
      if ( remove_expansion_rate_tensor() )
      {
        // txx
        set_val( txx_blk_idx(i,j,k), u_blk_idx(i,j,k),   -4.0/3.0);
        set_val( txx_blk_idx(i,j,k), u_blk_idx(i+1,j,k),  4.0/3.0);
        set_val( txx_blk_idx(i,j,k), v_blk_idx(i,j,k),    2.0/3.0);
        set_val( txx_blk_idx(i,j,k), v_blk_idx(i,j+1,k), -2.0/3.0);
        set_val( txx_blk_idx(i,j,k), w_blk_idx(i,j,k),    2.0/3.0);
        set_val( txx_blk_idx(i,j,k), w_blk_idx(i,j,k+1), -2.0/3.0);

        // tyy
        set_val( tyy_blk_idx(i,j,k), v_blk_idx(i,j,k),   -4.0/3.0);
        set_val( tyy_blk_idx(i,j,k), v_blk_idx(i,j+1,k),  4.0/3.0);
        set_val( tyy_blk_idx(i,j,k), u_blk_idx(i,j,k),    2.0/3.0);
        set_val( tyy_blk_idx(i,j,k), u_blk_idx(i+1,j,k), -2.0/3.0);
        set_val( tyy_blk_idx(i,j,k), w_blk_idx(i,j,k),    2.0/3.0);
        set_val( tyy_blk_idx(i,j,k), w_blk_idx(i,j,k+1), -2.0/3.0);

        // tzz
        set_val( tzz_blk_idx(i,j,k), w_blk_idx(i,j,k),   -4.0/3.0);
        set_val( tzz_blk_idx(i,j,k), w_blk_idx(i,j,k+1),  4.0/3.0);
        set_val( tzz_blk_idx(i,j,k), u_blk_idx(i,j,k),    2.0/3.0);
        set_val( tzz_blk_idx(i,j,k), u_blk_idx(i+1,j,k), -2.0/3.0);
        set_val( tzz_blk_idx(i,j,k), v_blk_idx(i,j,k),    2.0/3.0);
        set_val( tzz_blk_idx(i,j,k), v_blk_idx(i,j+1,k), -2.0/3.0);
      }
      else
      {
        // txx
        set_val( txx_blk_idx(i,j,k), u_blk_idx(i,j,k),   -2);
        set_val( txx_blk_idx(i,j,k), u_blk_idx(i+1,j,k),  2);

        // tyy
        set_val( tyy_blk_idx(i,j,k), v_blk_idx(i,j,k),   -2);
        set_val( tyy_blk_idx(i,j,k), v_blk_idx(i,j+1,k),  2);

        // tzz
        set_val( tzz_blk_idx(i,j,k), w_blk_idx(i,j,k),   -2);
        set_val( tzz_blk_idx(i,j,k), w_blk_idx(i,j,k+1),  2);
      }
    }
  }

  vit.setConstArray(myTyzIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if ( !isInSystem(vit.getValue()) )
      continue;

    // tyz
    set_val( tyz_blk_idx(i,j,k), v_blk_idx(i,j,k-1), -1);
    set_val( tyz_blk_idx(i,j,k), v_blk_idx(i,j,k),    1);
    set_val( tyz_blk_idx(i,j,k), w_blk_idx(i,j-1,k), -1);
    set_val( tyz_blk_idx(i,j,k), w_blk_idx(i,j,k),    1);
  }

  vit.setConstArray(myTxzIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if ( !isInSystem(vit.getValue()) )
      continue;

    // txz
    set_val( txz_blk_idx(i,j,k), u_blk_idx(i,j,k-1), -1);
    set_val( txz_blk_idx(i,j,k), u_blk_idx(i,j,k),    1);
    set_val( txz_blk_idx(i,j,k), w_blk_idx(i-1,j,k), -1);
    set_val( txz_blk_idx(i,j,k), w_blk_idx(i,j,k),    1);
  }

  vit.setConstArray(myTxyIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if ( !isInSystem(vit.getValue()) )
      continue;

    // txy
    set_val( txy_blk_idx(i,j,k), u_blk_idx(i,j-1,k), -1);
    set_val( txy_blk_idx(i,j,k), u_blk_idx(i,j,k),    1);
    set_val( txy_blk_idx(i,j,k), v_blk_idx(i-1,j,k), -1);
    set_val( txy_blk_idx(i,j,k), v_blk_idx(i,j,k),    1);
  }

  D.setFromTriplets( triplets.begin(), triplets.end() );

//  D *= 0.5;
}


template<typename T>
void
sim_stokesSolver<T>::buildGradientOperator(BlockMatrixType &G) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();
  auto set_val = [&](int row, int col, T val)
  {
    if ( isInSystem(row) )
      triplets.emplace_back( row, col, val );
  };

  // build gradient opperator on a per column basis
  UT_VoxelArrayIteratorI vit;
  vit.setConstArray(myCentralIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if ( !isInSystem(vit.getValue()) )
      continue;

    set_val( u_blk_idx(i,j,k),   p_blk_idx(i,j,k),  1);
    set_val( u_blk_idx(i+1,j,k), p_blk_idx(i,j,k), -1);
    set_val( v_blk_idx(i,j,k),   p_blk_idx(i,j,k),  1);
    set_val( v_blk_idx(i,j+1,k), p_blk_idx(i,j,k), -1);
    set_val( w_blk_idx(i,j,k),   p_blk_idx(i,j,k),  1);
    set_val( w_blk_idx(i,j,k+1), p_blk_idx(i,j,k), -1);
  }

  G.setFromTriplets( triplets.begin(), triplets.end() );
}

template<typename T>
void
sim_stokesSolver<T>::buildGhostFluidMatrix(
    const UT_VoxelArrayF &u_weights,
    const UT_VoxelArrayF &v_weights,
    const UT_VoxelArrayF &w_weights,
    BlockMatrixType &GF) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();

  auto add_gf_val = [&](int i, int j, int k, int in, int jn, int kn, int idx, T theta)
  {
    if ( p_blk_idx(in,jn,kn) == AIR || p_blk_idx(i,j,k) == AIR )
    {
      assert( theta );
      auto gf = (1.0-theta) / theta; // theta is guaranteed to be non-zero
      triplets.emplace_back( idx, idx, gf );
    }
  };

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&u_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = u_blk_idx(i,j,k);
    if (!isInSystem(idx) || vit.getValue() == 0.0f)
      continue;
    triplets.emplace_back(idx, idx, 1); // this matrix is the identity except for ghost fluid contributions
    add_gf_val(i,j,k, i-1,j,k, idx, u_weights(i,j,k));
  }
  vit.setConstArray(&v_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = v_blk_idx(i,j,k);
    if (!isInSystem(idx) || vit.getValue() == 0.0f)
      continue;
    triplets.emplace_back(idx, idx, 1);
    add_gf_val(i,j,k, i,j-1,k, idx, v_weights(i,j,k));
  }
  vit.setConstArray(&w_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = w_blk_idx(i,j,k);
    if (!isInSystem(idx) || vit.getValue() == 0.0f)
      continue;
    triplets.emplace_back(idx, idx, 1);
    add_gf_val(i,j,k, i,j,k-1, idx, w_weights(i,j,k));
  }

  GF.setFromTriplets( triplets.begin(), triplets.end() );
}

template<typename T>
void
sim_stokesSolver<T>::buildSumNeighboursOperator(BlockMatrixType &N) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();
  auto set_val = [&](int row, int col, T val)
  {
    if ( isInSystem(col) )
      triplets.emplace_back( row, col, val );
  };

  // build gradient opperator on a per column basis
  UT_VoxelArrayIteratorI vit;
  vit.setConstArray(myCentralIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if ( !isInSystem(vit.getValue()) )
      continue;

    set_val( p_blk_idx(i,j,k), u_blk_idx(i,j,k),   1);
    set_val( p_blk_idx(i,j,k), u_blk_idx(i+1,j,k), 1);
    set_val( p_blk_idx(i,j,k), v_blk_idx(i,j,k),   1);
    set_val( p_blk_idx(i,j,k), v_blk_idx(i,j+1,k), 1);
    set_val( p_blk_idx(i,j,k), w_blk_idx(i,j,k),   1);
    set_val( p_blk_idx(i,j,k), w_blk_idx(i,j,k+1), 1);
  }

  N.setFromTriplets( triplets.begin(), triplets.end() );
}

template<typename T>
template<bool INVERSE>
void
sim_stokesSolver<T>::buildViscosityMatrix(
    const SIM_RawField &viscfield, BlockMatrixType &M) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  const UT_VoxelArrayF &viscosity = *viscfield.field();

  triplets.clear();

  UT_VoxelArrayIteratorI vit;
  vit.setConstArray(myCentralIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if ( !isInSystem(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    float visc = viscosity(i,j,k);
    if ( visc > 0 )
    {
      if ( reduced_stress_tensor() )
      {bool
        T diag, offdiag;
        if ( INVERSE )
        {
          diag = 1.0/visc;
          offdiag = 0.5/visc;
        }
        else
        {
          diag = visc*4.0/3.0;
          offdiag = -visc*2.0/3.0;
        }

        // txx
        triplets.emplace_back(txx_blk_idx(i,j,k), txx_blk_idx(i,j,k), diag);
        triplets.emplace_back(txx_blk_idx(i,j,k), tyy_blk_idx(i,j,k), offdiag);

        // tyy
        triplets.emplace_back(tyy_blk_idx(i,j,k), tyy_blk_idx(i,j,k), diag);
        triplets.emplace_back(tyy_blk_idx(i,j,k), txx_blk_idx(i,j,k), offdiag);
      }
      else
      {
        T val = INVERSE ? 1.0/visc : visc;
        // txx, tyy, tzz
        triplets.emplace_back(txx_blk_idx(i,j,k), txx_blk_idx(i,j,k), val);
        triplets.emplace_back(tyy_blk_idx(i,j,k), tyy_blk_idx(i,j,k), val);
        triplets.emplace_back(tzz_blk_idx(i,j,k), tzz_blk_idx(i,j,k), val);
      }
    }
    else
      assert("bad viscosity");
  }

  UT_VoxelProbeAverage<float, 0, -1, -1> probe_yz;
  probe_yz.setArray(&viscosity);
  vit.setConstArray(myTyzIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if( !isInSystem(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    // tyz
    probe_yz.setIndex(vit);
    float visc = probe_yz.getValue();
    //assert( viscosity.lerpVoxel(i,j,k,0,-0.5,-0.5) == visc );
    if ( visc > 0 )
      triplets.emplace_back(tyz_blk_idx(i,j,k), tyz_blk_idx(i,j,k),  INVERSE ? 1.0/visc : visc);
    else
      assert("bad viscosity");
  }

  UT_VoxelProbeAverage<float, -1, 0, -1> probe_xz;
  probe_xz.setArray(&viscosity);
  vit.setConstArray(myTxzIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if ( !isInSystem(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    // txz
    probe_xz.setIndex(vit);
    float visc = probe_xz.getValue();
    //assert( viscosity.lerpVoxel(i,j,k,-0.5,0,-0.5) == visc );
    if ( visc > 0 )
      triplets.emplace_back(txz_blk_idx(i,j,k), txz_blk_idx(i,j,k), INVERSE ? 1.0/visc : visc);
    else
      assert("bad viscosity");
  }

  UT_VoxelProbeAverage<float, -1, -1, 0> probe_xy;
  probe_xy.setArray(&viscosity);
  vit.setConstArray(myTxyIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if ( !isInSystem(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    // txy
    probe_xy.setIndex(vit);
    float visc = probe_xy.getValue();
    //assert( viscosity.lerpVoxel(i,j,k,-0.5,-0.5,0) == visc );
    if ( visc > 0 )
      triplets.emplace_back(txy_blk_idx(i,j,k), txy_blk_idx(i,j,k), INVERSE ? 1.0/visc : visc);
    else
      assert("bad viscosity");
  }

  M.setFromTriplets( triplets.begin(), triplets.end() );
}

template<typename T>
void
sim_stokesSolver<T>::buildDensityMatrix(
    const SIM_RawField &densfield,
    BlockMatrixType &P) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  const UT_VoxelArrayF &density = *densfield.field();

  // NOTE: if inverse is required, consider what happens at thet solid boundary.
  // Currently it makes sense for invalid (collision) cells to have infinite
  // density, and thus 1.0/density -> 0, so it's ok to skip those.
  triplets.clear();

  UT_VoxelArrayIteratorI vit;

  UT_VoxelProbeAverage<float,-1,0,0> probe_x;
  probe_x.setArray(&density);
  vit.setConstArray(myUIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if ( !isInSystem(vit.getValue()) || isCollision(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    probe_x.setIndex(vit);
    float dens = probe_x.getValue();
    //assert( density.lerpVoxel(i,j,k,-0.5,0,0) == dens )
    if ( dens > 0 )
      triplets.emplace_back( u_blk_idx(i,j,k), u_blk_idx(i,j,k),  1.0/dens);
    else
      assert("bad density");
  }

  UT_VoxelProbeAverage<float,0,-1,0> probe_y;
  probe_y.setArray(&density);
  vit.setConstArray(myVIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if( !isInSystem(vit.getValue()) || isCollision(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    probe_y.setIndex(vit);
    float dens = probe_y.getValue();
    //assert( density.lerpVoxel(i,j,k,0,-0.5,0) != dens )
    if ( dens > 0 )
      triplets.emplace_back(v_blk_idx(i,j,k), v_blk_idx(i,j,k),  1.0/dens);
    else
      assert("bad density");
  }

  UT_VoxelProbeAverage<float,0,0,-1> probe_z;
  probe_z.setArray(&density);
  vit.setConstArray(myWIndex.field());
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    if( !isInSystem(vit.getValue()) || isCollision(vit.getValue()) )
      continue;

    int i = vit.x(), j = vit.y(), k = vit.z();

    probe_z.setIndex(vit);
    float dens = probe_z.getValue();
    //assert( density.lerpVoxel(i,j,k,0,0,-0.5) != dens )
    if ( dens > 0 )
      triplets.emplace_back(w_blk_idx(i,j,k), w_blk_idx(i,j,k),  1.0/dens);
    else
      assert("bad density");
  }

  P.setFromTriplets( triplets.begin(), triplets.end() );
}

template<typename T>
void
sim_stokesSolver<T>::buildPressureWeightMatrix(
    const UT_VoxelArrayF &c_weights,
    BlockMatrixType& W) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&c_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = p_blk_idx(i,j,k);
    if (!isInSystem(idx) || vit.getValue() == 0.0f)
      continue;
    triplets.emplace_back(idx, idx, vit.getValue());
  }

  W.setFromTriplets(triplets.begin(), triplets.end());
}

template<typename T>
template<bool INVERSE>
void
sim_stokesSolver<T>::buildVelocityWeightMatrix(
    const UT_VoxelArrayF &u_weights,
    const UT_VoxelArrayF &v_weights,
    const UT_VoxelArrayF &w_weights,
    BlockMatrixType& W) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&u_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(u_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    auto idx = u_blk_idx(i,j,k);
    triplets.emplace_back(idx, idx, INVERSE ? 1.0/vit.getValue() : vit.getValue());
  }
  vit.setConstArray(&v_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(v_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    auto idx = v_blk_idx(i,j,k);
    triplets.emplace_back(idx, idx, INVERSE ? 1.0/vit.getValue() : vit.getValue());
  }
  vit.setConstArray(&w_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(w_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    auto idx = w_blk_idx(i,j,k);
    triplets.emplace_back(idx, idx, INVERSE ? 1.0/vit.getValue() : vit.getValue());
  }

  W.setFromTriplets(triplets.begin(), triplets.end());
}

template<typename T>
template<bool INVERSE>
void
sim_stokesSolver<T>::buildStressWeightMatrix(
    const UT_VoxelArrayF &c_weights,
    const UT_VoxelArrayF &ex_weights,
    const UT_VoxelArrayF &ey_weights,
    const UT_VoxelArrayF &ez_weights,
    BlockMatrixType &W) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  triplets.clear();
  auto set_val = [&](int idx, T val)
  {
      triplets.emplace_back(idx, idx, val);
  };

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&c_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(p_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    auto val = INVERSE ? 1.0/vit.getValue() : vit.getValue();
    set_val(txx_blk_idx(i,j,k), val);
    set_val(tyy_blk_idx(i,j,k), val);
    if (!reduced_stress_tensor())
      set_val(tzz_blk_idx(i,j,k), val);
  }
  vit.setConstArray(&ex_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(tyz_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    set_val(tyz_blk_idx(i,j,k), INVERSE ? 1.0/vit.getValue() : vit.getValue());
  }
  vit.setConstArray(&ey_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(txz_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    set_val(txz_blk_idx(i,j,k), INVERSE ? 1.0/vit.getValue() : vit.getValue());
  }
  vit.setConstArray(&ez_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(txy_blk_idx(i,j,k)) || vit.getValue() == 0.0f)
      continue;
    set_val(txy_blk_idx(i,j,k), INVERSE ? 1.0/vit.getValue() : vit.getValue());
  }

  W.setFromTriplets(triplets.begin(), triplets.end());
}



// PRE:  input is the surface tension pressure field
// POST: output is the vector aligned with velocity values at cell interfaces
// containing the ghost pressures where one neighbour is inside and one outside
// of the liquid. Values where there is no ghost pressure are zero.
template<typename T>
auto
sim_stokesSolver<T>::buildGhostFluidSurfaceTensionPressureVector(
    const SIM_RawField * const* surf_weights,
    const SIM_RawField & surfp) const -> BlockVectorType{
  const UT_VoxelArrayF &u_weights = *surf_weights[4]->field();
  const UT_VoxelArrayF &v_weights = *surf_weights[5]->field();
  const UT_VoxelArrayF &w_weights = *surf_weights[6]->field();
  const UT_VoxelArrayF &sp = *(surfp.field());

  BlockVectorType ust(myNumVelocityVars);
  ust.setZero();

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&u_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto uidx = u_blk_idx(i,j,k);
    if (!isInSystem(uidx) || vit.getValue() == 0.0f)
      continue;
    ust[uidx] = ghostFluidSurfaceTensionPressure<0>(i,j,k, vit.getValue(), sp);
  }
  vit.setConstArray(&v_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto vidx = v_blk_idx(i,j,k);
    if (!isInSystem(vidx) || vit.getValue() == 0.0f)
      continue;
    ust[vidx] = ghostFluidSurfaceTensionPressure<1>(i,j,k, vit.getValue(), sp);
  }
  vit.setConstArray(&w_weights);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto widx = w_blk_idx(i,j,k);
    if (!isInSystem(widx) || vit.getValue() == 0.0f)
      continue;
    ust[widx] = ghostFluidSurfaceTensionPressure<2>(i,j,k, vit.getValue(), sp);
  }

  return ust;
}

template<typename T>
auto
sim_stokesSolver<T>::buildSolidVelocityVector(
    const SIM_RawField * const *vel) const -> BlockVectorType{
  BlockVectorType uout(myNumVelocityVars);
  uout.setZero();

    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  const UT_VoxelArrayF &u = *(vel[0]->field());
  const UT_VoxelArrayF &v = *(vel[1]->field());
  const UT_VoxelArrayF &w = *(vel[2]->field());

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&u);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = u_blk_idx(i,j,k);
    if (!isCollision(idx))
      continue;
    uout[idx] = vit.getValue();
  }
  vit.setConstArray(&v);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = v_blk_idx(i,j,k);
    if (!isCollision(idx))
      continue;
    uout[idx] = vit.getValue();
  }
  vit.setConstArray(&w);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = w_blk_idx(i,j,k);
    if (!isCollision(idx))
      continue;
    uout[idx] = vit.getValue();
  }
  return uout;
}

template<typename T>
auto
sim_stokesSolver<T>::buildSurfaceTensionPressureVector(const SIM_RawField &surfp) const -> BlockVectorType{
  BlockVectorType pout(myNumPressureVars);
  pout.setZero();

    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  const UT_VoxelArrayF &sp = *(surfp.field());

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&sp);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    if (!isInSystem(p_blk_idx(i,j,k)))
      continue;
    pout[p_blk_idx(i,j,k)] = vit.getValue();
  }
  return pout;
}

template<typename T>
auto
sim_stokesSolver<T>::buildVelocityVector(
    const SIM_VectorField &vel,
    const SIM_RawField * const * colvel) const -> BlockVectorType{
  BlockVectorType ustar(myNumVelocityVars);
  ustar.setZero();

    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  const UT_VoxelArrayF &u_col = *colvel[0]->field();
  const UT_VoxelArrayF &v_col = *colvel[1]->field();
  const UT_VoxelArrayF &w_col = *colvel[2]->field();

  const UT_VoxelArrayF &u = *vel.getField(0)->field();
  const UT_VoxelArrayF &v = *vel.getField(1)->field();
  const UT_VoxelArrayF &w = *vel.getField(2)->field();

  UT_VoxelArrayIteratorF vit;
  vit.setConstArray(&u);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = u_blk_idx(i,j,k);
    if (!isInSystem(idx))
      continue;
    ustar[idx] = isCollision(idx) ? u_col(i,j,k) : vit.getValue();
  }
  vit.setConstArray(&v);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = v_blk_idx(i,j,k);
    if (!isInSystem(idx))
      continue;
    ustar[idx] = isCollision(idx) ? v_col(i,j,k) : vit.getValue();
  }
  vit.setConstArray(&w);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    auto idx = w_blk_idx(i,j,k);
    if (!isInSystem(idx))
      continue;
    ustar[idx] = isCollision(idx) ? w_col(i,j,k) : vit.getValue();
  }
  return ustar;
}

// returns false if nans or infinite values are found in the given matrix
template<typename MatrixType>
bool
isMatrixValid(const MatrixType& A){
  for ( int k = 0; k < A.outerSize(); ++k )
  {
    for ( typename MatrixType::InnerIterator it(A, k); it; ++it )
    {
      if ( SYSisNan(it.value()) || !SYSisFinite(it.value()) )
      {
          std::cerr << "(" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;

          return false;
      }
    }
  }
  return true;
}

// returns false if symmetric matrix A contains any rows (or columns) of all zeros
template<typename MatrixType>
bool
isMatrixPruned(const MatrixType& A){
  for (int k = 0; k < A.outerSize(); ++k )
  {
    typename MatrixType::InnerIterator it(A,k);
    if (!it)
    {
      return false;
    }
  }
  return true;
}

template<typename T>
void sim_stokesSolver<T>::assembleBlockSystem(
    const BlockMatrixType& WLp,
    const BlockMatrixType& WLuinv,
    const BlockMatrixType& WFu,
    const BlockMatrixType& WLt,
    const BlockMatrixType& WFt,
    const BlockMatrixType& G,
    const BlockMatrixType& D,
    const BlockMatrixType& Pinv,
    const BlockMatrixType& Minv,
    BlockMatrixType& Ap,
    BlockMatrixType& Bp,
    BlockMatrixType& Hp,
    BlockMatrixType& At,
    BlockMatrixType& Bt,
    BlockMatrixType& Ht) const{
  // Original expressions
  //Ap = dt*WLp*G.transpose()*Pinv*WLuinv*WFu*G*WLp;
  //Bp = dx*WLp*G.transpose()*WFu;  // B where b = B*ustar
  //Hp = (dt/dx)*WLuinv*Pinv*G*WLp;  // H update matrix
  //At = (dx*dx*0.5)*Minv*WLt*WFt + dt*WLt*D*Pinv*WLuinv*WFu*D.transpose()*WLt;
  //Bt = dx*WLt*D*WFu;  // B where b = B*ustar
  //Ht = (dt/dx)*WLuinv*Pinv*D.transpose()*WLt;  // H update matrix

  // Optimized matrix operations
  BlockMatrixType PinvDtransposeWLt, WLtDWFu, WLtD, GWLp, PinvGWLp, WLpGtransposeWFu;
  WLtD = WLt*D;
  WLtDWFu = WLtD*WFu;
  PinvDtransposeWLt = Pinv*WLtD.transpose();
  GWLp = G*WLp;
  WLpGtransposeWFu = GWLp.transpose()*WFu;
  PinvGWLp = Pinv*GWLp;
  At = (dx*dx)*Minv*WLt*WFt + dt*WLtDWFu*WLuinv*PinvDtransposeWLt;
  Bt = dx*WLtDWFu;  // B where b = B*ustar
  Ht = (dt/dx)*WLuinv*PinvDtransposeWLt;  // H update matrix

  Ap = dt*WLpGtransposeWFu*WLuinv*PinvGWLp;
  Bp = dx*WLpGtransposeWFu;     // B where b = B*ustar
  Hp = (dt/dx)*WLuinv*PinvGWLp; // H update matrix
}

template<typename T>
void sim_stokesSolver<T>::assembleStressVelocitySystem(
    const BlockMatrixType& WLt,
    const BlockMatrixType& WLu,
    const BlockMatrixType& WFtinv,
    const BlockMatrixType& WFu,
    const BlockMatrixType& D,
    const BlockMatrixType& Pinv,
    const BlockMatrixType& M,
    BlockMatrixType& A,
    BlockMatrixType& B) const{
  auto dx2 = dx*dx;
  A = dx2*WLu*WFu + 0.5*dt*Pinv*WFu*D.transpose()*M*WFtinv*WLt*D*WFu;
  B = dx2*WLu*WFu;  // B where b = B*ustar
}

template<typename T>
void sim_stokesSolver<T>::buildSystemBlockwise(
    BlockMatrixType &matrix,
    BlockVectorType &rhs,
    BlockMatrixType &H,
    BlockVectorType &ust,
    const BlockVectorType &uold,
    const SIM_RawField & surf,
    const SIM_RawField * const* surf_weights,
    const SIM_RawField * const* col_weights,
    const SIM_RawField &viscfield,
    const SIM_RawField &densfield,
    const SIM_RawField * const* solid_vel,
    const SIM_RawField & surf_pres) const{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif
  const UT_VoxelArrayF &c_vol_liquid = *surf_weights[0]->field();
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

  // here we will take another approach and construct the system blockwize:
  // building Att, Atp and App separately.

  BlockMatrixType D(myNumStressVars, myNumVelocityVars);
  BlockMatrixType G(myNumVelocityVars, myNumPressureVars);
  BlockMatrixType Minv(myNumStressVars, myNumStressVars);
  BlockMatrixType Pinv(myNumVelocityVars, myNumVelocityVars);

  BlockMatrixType WLp(myNumPressureVars, myNumPressureVars);
  BlockMatrixType WFp(myNumPressureVars, myNumPressureVars);
  BlockMatrixType WLt(myNumStressVars, myNumStressVars);
  BlockMatrixType WFt(myNumStressVars, myNumStressVars);
  BlockMatrixType WLuinv(myNumVelocityVars, myNumVelocityVars);
  BlockMatrixType WFu(myNumVelocityVars, myNumVelocityVars);

  buildDeformationRateOperator(D);
  buildGradientOperator(G);
  buildViscosityMatrix<true>(viscfield, Minv);
  buildDensityMatrix(densfield, Pinv);
  buildPressureWeightMatrix(c_vol_liquid, WLp);

  buildStressWeightMatrix<false>(c_vol_liquid, ex_vol_liquid, ey_vol_liquid, ez_vol_liquid, WLt);
  buildVelocityWeightMatrix<true>(u_vol_liquid, v_vol_liquid, w_vol_liquid, WLuinv);
  buildStressWeightMatrix<false>(c_vol_fluid, ex_vol_fluid, ey_vol_fluid, ez_vol_fluid, WFt);
  buildVelocityWeightMatrix<false>(u_vol_fluid, v_vol_fluid, w_vol_fluid, WFu);
  buildPressureWeightMatrix(c_vol_fluid, WFp);

  BlockMatrixType App;// = dt*WLp*G.transpose()*Pinv*WLuinv*WFu*G*WLp;
  BlockMatrixType Att;// = (dx*dx*0.5)*Minv*WLt*WFt + dt*WLt*D*Pinv*WLuinv*WFu*D.transpose()*WLt;
  BlockMatrixType Bp, Bt, Ht, Hp, Atp;

  assembleBlockSystem(WLp, WLuinv, WFu, WLt, WFt, G, D, Pinv, Minv, App, Bp, Hp, Att, Bt, Ht);
  Atp = dt*WLt*D*Pinv*WLuinv*WFu*G*WLp;
  // Print matrices in dense form
  //MatrixX<T> Attdense, Appdense, Atpdense;
  //Attdense = MatrixX<T>(Att);
  //Appdense = MatrixX<T>(App);
  //Atpdense = MatrixX<T>(Atp);
  //std::cerr << "Att = " << Attdense << std::endl;
  //std::cerr << "App = " << Appdense << std::endl;
  //std::cerr << "Atp = " << Atpdense << std::endl;
  //
  //std::cerr << "G^T = \n" << G.transpose() << std::endl;
  //BlockMatrixType GtG = G.transpose() * G;
  //std::cerr << "G^T*G = \n" << GtG << std::endl;
  //std::cerr << "App = \n" << App << std::endl;
  //std::cerr << "Pinv = \n" << Pinv << std::endl;
  //std::cerr << "WLp = \n" << WLp << std::endl;

  assert(isMatrixValid(App));
  assert(isMatrixValid(Att));
  assert(isMatrixValid(Atp));

  auto elts = getNumStokesVars();

  // moving boundary term
  BlockVectorType ubc = buildSolidVelocityVector(solid_vel);

    //////
    // TODO: verify that the ust from the two lines below is the same as the ust
    // we actually use. To do this we have to modify
    // buildSurfaceTensionPressureVector to include pressures at air cells and
    // exclude other interior pressures (not near the surface). This will show that
    // we can compute ghost pressures (gfst) using the formula:
    //    (WLu*G - G*WLp)*pbc
    // in the same fashion we compute the boundary velocity condition
    // (We can also do this on paper, what what's the fun in that :P)
    //////
    //  BlockVectorType pbc = buildSurfaceTensionPressureVector(surf_pres);
    //  ust = buildSurfaceTensionRHSAlt(surf_weights, densfield, pbc);
    //////

  // surface tension term
  BlockVectorType gfst = buildGhostFluidSurfaceTensionPressureVector(surf_weights, surf_pres);
  ust = buildSurfaceTensionRHS(surf_weights, densfield, gfst);
  //ust = dt*Pinv*(G - WLuinv*G*WLp)*pbc;

  rhs.resize(elts);
  rhs.setZero();
  rhs << Bp*uold - dx*WLp*(G.transpose()*WFu - WFp*G.transpose())*ubc + WLp*G.transpose()*WFu*ust,
         Bt*uold - dx*WLt*(D*WFu - WFt*D)*ubc + WLt*D*WFu*ust;

  triplets.clear();
  matrix.resize(elts,elts);

  // Copy App matrix
  for ( int k = 0; k < App.outerSize(); ++k )
    for ( typename BlockMatrixType::InnerIterator it(App,k); it; ++it )
      triplets.emplace_back(it.row(), it.col(), it.value());

  // Copy Att matrix
  for ( int k = 0; k < Att.outerSize(); ++k )
    for ( typename BlockMatrixType::InnerIterator it(Att,k); it; ++it )
    {
      int row = it.row() + myNumPressureVars;
      int col = it.col() + myNumPressureVars;
      triplets.emplace_back(row, col, it.value());
    }

  // Copy Atp matrix
  for ( int k = 0; k < Atp.outerSize(); ++k )
    for ( typename BlockMatrixType::InnerIterator it(Atp,k); it; ++it )
    {
      int row = it.row() + myNumPressureVars;
      int col = it.col();
      triplets.emplace_back(row, col, it.value());
      triplets.emplace_back(col, row, it.value()); // transpose
    }

  matrix.setFromTriplets(triplets.begin(), triplets.end());

  triplets.clear();

  // Copy Hp and Ht matrices matrix
  for ( int k = 0; k < Hp.outerSize(); ++k )
    for ( typename BlockMatrixType::InnerIterator it(Hp,k); it; ++it )
      triplets.emplace_back(it.row(), it.col(), it.value());

  for ( int k = 0; k < Ht.outerSize(); ++k )
    for ( typename BlockMatrixType::InnerIterator it(Ht,k); it; ++it )
      triplets.emplace_back(it.row(), it.col()+myNumPressureVars, it.value());

  H.resize(myNumVelocityVars, myNumPressureVars+myNumStressVars);
  H.setFromTriplets(triplets.begin(), triplets.end());
}

template<typename T>
auto
sim_stokesSolver<T>::buildSurfaceTensionRHSAlt(
    const SIM_RawField * const* surf_weights,
    const SIM_RawField & densfield,
    const BlockVectorType &pbc) const -> BlockVectorType{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif

  const UT_VoxelArrayF &c_vol_liquid = *surf_weights[0]->field();
  const UT_VoxelArrayF &u_vol_liquid = *surf_weights[4]->field();
  const UT_VoxelArrayF &v_vol_liquid = *surf_weights[5]->field();
  const UT_VoxelArrayF &w_vol_liquid = *surf_weights[6]->field();

  BlockMatrixType G(myNumVelocityVars, myNumPressureVars);
  BlockMatrixType Pinv(myNumVelocityVars, myNumVelocityVars);
  BlockMatrixType WLp(myNumPressureVars, myNumPressureVars);
  BlockMatrixType WLuinv(myNumVelocityVars, myNumVelocityVars);

  buildGradientOperator(G);
  buildDensityMatrix(densfield, Pinv);
  buildVelocityWeightMatrix<true>(u_vol_liquid, v_vol_liquid, w_vol_liquid, WLuinv);
  buildPressureWeightMatrix(c_vol_liquid, WLp);

  return dt*Pinv*(G - WLuinv*G*WLp)*pbc;
}

template<typename T>
auto
sim_stokesSolver<T>::buildSurfaceTensionRHS(
    const SIM_RawField * const* surf_weights,
    const SIM_RawField & densfield,
    const BlockVectorType &ust) const -> BlockVectorType{
    #ifndef BLOCKWISE_STOKES
    assert(myScheme != STOKES);
    #endif

  const UT_VoxelArrayF &u_vol_liquid = *surf_weights[4]->field();
  const UT_VoxelArrayF &v_vol_liquid = *surf_weights[5]->field();
  const UT_VoxelArrayF &w_vol_liquid = *surf_weights[6]->field();

  BlockMatrixType Pinv(myNumVelocityVars, myNumVelocityVars);
  BlockMatrixType WLuinv(myNumVelocityVars, myNumVelocityVars);

  buildDensityMatrix(densfield, Pinv);
  buildVelocityWeightMatrix<true>(u_vol_liquid, v_vol_liquid, w_vol_liquid, WLuinv);

  // dt*Pinv*(G - WLuinv*G*WLp)*pbc;
  return (-dt) * Pinv * WLuinv * ust;
}

template<typename T>
void sim_stokesSolver<T>::buildDecoupledSystem(
    BlockMatrixType &At, BlockMatrixType &Bt, BlockMatrixType &Ht,
    BlockMatrixType &Ap, BlockMatrixType &Bp, BlockMatrixType &Hp,
    const SIM_RawField & surf,
    const SIM_RawField * const* surf_weights,
    const SIM_RawField * const* col_weights,
    const SIM_RawField &viscfield,
    const SIM_RawField &densfield,
    const SIM_RawField * const* solid_vel,
    const SIM_RawField & surf_pres) const{
  assert(myScheme != STOKES);
  const UT_VoxelArrayF &c_vol_liquid = *surf_weights[0]->field();
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

  BlockMatrixType D(myNumStressVars, myNumVelocityVars);
  BlockMatrixType G(myNumVelocityVars, myNumPressureVars);
  BlockMatrixType Minv(myNumStressVars, myNumStressVars);
  BlockMatrixType Pinv(myNumVelocityVars, myNumVelocityVars);

  BlockMatrixType WLt(myNumStressVars, myNumStressVars);
  BlockMatrixType WFt(myNumStressVars, myNumStressVars);
  BlockMatrixType WLp(myNumPressureVars, myNumPressureVars);
  BlockMatrixType WLuinv(myNumVelocityVars, myNumVelocityVars);
  BlockMatrixType WFu(myNumVelocityVars, myNumVelocityVars);

  buildDeformationRateOperator(D);
  buildGradientOperator(G);
  buildViscosityMatrix<true>(viscfield, Minv);
  buildDensityMatrix(densfield, Pinv);

  buildStressWeightMatrix<false>(c_vol_liquid, ex_vol_liquid, ey_vol_liquid, ez_vol_liquid, WLt);
  buildPressureWeightMatrix(c_vol_liquid, WLp);
  buildVelocityWeightMatrix<true>(u_vol_liquid, v_vol_liquid, w_vol_liquid, WLuinv);
  buildStressWeightMatrix<false>(c_vol_fluid, ex_vol_fluid, ey_vol_fluid, ez_vol_fluid, WFt);
  buildVelocityWeightMatrix<false>(u_vol_fluid, v_vol_fluid, w_vol_fluid, WFu);

  assembleBlockSystem(WLp, WLuinv, WFu, WLt, WFt, G, D, Pinv, Minv, Ap, Bp, Hp, At, Bt, Ht);
}

template<typename T>
void sim_stokesSolver<T>::buildPressureOnlySystem(
    BlockMatrixType &A, BlockMatrixType &B, BlockMatrixType &H,
    const SIM_RawField & surf,
    const SIM_RawField * const* surf_weights,
    const SIM_RawField * const* col_weights,
    const SIM_RawField &densfield,
    const SIM_RawField * const* solid_vel,
    const SIM_RawField & surf_pres) const{
  assert(myScheme != STOKES);

  const UT_VoxelArrayF &u_vol_liquid = *surf_weights[4]->field();
  const UT_VoxelArrayF &v_vol_liquid = *surf_weights[5]->field();
  const UT_VoxelArrayF &w_vol_liquid = *surf_weights[6]->field();

  const UT_VoxelArrayF &u_vol_fluid = *col_weights[4]->field();
  const UT_VoxelArrayF &v_vol_fluid = *col_weights[5]->field();
  const UT_VoxelArrayF &w_vol_fluid = *col_weights[6]->field();

  BlockMatrixType G(myNumVelocityVars, myNumPressureVars);
  BlockMatrixType Pinv(myNumVelocityVars, myNumVelocityVars);

  BlockMatrixType WLuinv(myNumVelocityVars, myNumVelocityVars);
  BlockMatrixType WFu(myNumVelocityVars, myNumVelocityVars);

  buildGradientOperator(G);

  buildDensityMatrix(densfield, Pinv);
  buildVelocityWeightMatrix<true>(u_vol_liquid, v_vol_liquid, w_vol_liquid, WLuinv);
  buildVelocityWeightMatrix<false>(u_vol_fluid, v_vol_fluid, w_vol_fluid, WFu);

    #if 0
    // TODO: figure out why GF is different than WLuinv. I'm not satisfied in
    // knowing that they are almost the same and WLuinv works.
    BlockMatrixType GF(myNumVelocityVars, myNumVelocityVars);
    buildGhostFluidMatrix(u_vol_liquid, v_vol_liquid, w_vol_liquid, GF);
    for ( int k = 0; k < myNumVelocityVars; ++k )
    {
        for (typename BlockMatrixType::InnerIterator it(GF,k); it; ++it)
        {
        auto diff = it.value() - WLuinv.coeff(it.row(), it.col());
        if ( 0 && diff )
        {
            std::cerr << "GF("<< it.row() << ", " << it.col() << " = "
            << it.value() << " vs. "
            << " WLuinv = " << WLuinv.coeff(it.row(), it.col()) << "; diff = " << diff << std::endl;
        }
        }
    }
    #endif

  BlockMatrixType GTWFu = G.transpose()*WFu;
  BlockMatrixType PinvWLuinvG = Pinv*WLuinv*G;
  A = dt*GTWFu*PinvWLuinvG;
  B = dx*GTWFu;             // B where b = B*ustar
  H = (dt/dx)*PinvWLuinvG;  // H update matrix
}

template<typename T>
void sim_stokesSolver<T>::buildViscositySystem(
    BlockMatrixType &Au, BlockMatrixType &Bu,
    const SIM_RawField * const* surf_weights,
    const SIM_RawField * const* col_weights,
    const SIM_RawField &viscfield,
    const SIM_RawField &densfield,
    const SIM_RawField * const* solid_vel) const{
  assert(myScheme != STOKES);
  const UT_VoxelArrayF &c_vol_liquid = *surf_weights[0]->field();
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

  BlockMatrixType D(myNumStressVars, myNumVelocityVars);
  BlockMatrixType G(myNumVelocityVars, myNumPressureVars);
  BlockMatrixType M(myNumStressVars, myNumStressVars);
  BlockMatrixType Pinv(myNumVelocityVars, myNumVelocityVars);

  BlockMatrixType WLt(myNumStressVars, myNumStressVars);
  BlockMatrixType WFtinv(myNumStressVars, myNumStressVars);
  BlockMatrixType WFu(myNumVelocityVars, myNumVelocityVars);
  BlockMatrixType WLu(myNumVelocityVars, myNumVelocityVars);

  buildDeformationRateOperator(D);
  buildViscosityMatrix<false>(viscfield, M);
  buildDensityMatrix(densfield, Pinv);

  buildStressWeightMatrix<false>(c_vol_liquid, ex_vol_liquid, ey_vol_liquid, ez_vol_liquid, WLt);
  buildVelocityWeightMatrix<false>(u_vol_liquid, v_vol_liquid, w_vol_liquid, WLu);
  buildStressWeightMatrix<true>(c_vol_fluid, ex_vol_fluid, ey_vol_fluid, ez_vol_fluid, WFtinv);
  buildVelocityWeightMatrix<false>(u_vol_fluid, v_vol_fluid, w_vol_fluid, WFu);

  assembleStressVelocitySystem(WLt, WLu, WFtinv, WFu, D, Pinv, M, Au, Bu);
}

template<typename T>
void
sim_stokesSolver<T>::removeNullSpace(const MatrixType &matrix, const VectorType &rhs) const
    {
    //For all-closed domains there may be a constant pressure nullspace
    //The solver can handle it if we suggest what the nullspace may be.
    /*
    std::vector<double> pressure_nullspace(rhs.size(),0);
    for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
        pressure_nullspace[p_idx(i,j,k)] = c_valid(i,j,k)?1:0;
    }

    */
    //  assert( matrix.m == matrix.n );

    //  Eigen::SparseMatrix<double> M(matrix.m, matrix.n);
    //  std::vector<Triplet<T>> triplets;
    //
    //  for(unsigned int i=0; i<matrix.m; ++i){
    //    for(unsigned int k=0; k<matrix.index[i].size(); ++k){
    //      triplets.push_back(Triplet<T>(i,matrix.index[i][k], matrix.value[i][k]));
    //    }
    //  }
    //  M.setFromTriplets(triplets.begin(), triplets.end());

    //for ( unsigned int i = 0; i < pressure_nullspace.size(); ++i )
    //{
    //  if ( !pressure_nullspace[i] )
    //    continue;

    //  rhs( i ) = 0.0;
    //  M.coeffRef(i,i) = 1.0;
    //  for ( Eigen::SparseMatrix<double>::InnerIterator it( M, i ); it; ++it )
    //  {
    //    auto& mtx_val = it.valueRef();
    //    if ( it.row() == it.col() )
    //      continue;
    //    if ( mtx_val == 0.0)
    //      continue;

    //    mtx_val = 0.0;
    //    if ( M.coeff( it.col(), it.row() ) != 0.0 )
    //        M.coeffRef( it.col(), it.row() ) = 0.0;
    //  }
    //}

    //std::cout << "system size: " << rhs.size() << "\n";
    //int sum =0;
    //for ( auto v : pressure_nullspace )
    //  sum += v;
    //std::cout << "nullspace size: " << sum << "\n";

    }

// remove zero rows and columns from the system
template<typename T>
void
sim_stokesSolver<T>::pruneSystem(
    // input
    const BlockMatrixType &A,
    const BlockVectorType &b,
    // output
    MatrixType            &newA, // colmajor
    VectorType            &newb,
    UT_ExintArray         &to_original) const{
  assert(!BlockMatrixType::IsRowMajor);
  to_original.clear();
  UT_ExintArray to_new(A.outerSize(), A.outerSize());
  to_new.constant(-1);

  for (int k = 0; k < A.outerSize(); ++k)
  {
    typename BlockMatrixType::InnerIterator it(A,k);
    if ( it )
    {
      to_new[k] = to_original.size();
      to_original.append(k);
    }
  }

  auto new_size = to_original.size();
  newb.init(0, new_size-1);
  newA.init(new_size,/* nonzeros = */29);
  std::vector<int> rowidx(new_size, 0); // per row
  for ( int k = 0; k < new_size; ++k )
  {
    auto orig_k = to_original[k];
    newb(k) = b[orig_k];
    // collect non-zeros for the new matrix
    for (typename BlockMatrixType::InnerIterator it(A,orig_k); it; ++it)
    {
      auto new_row = to_new[it.row()];
      assert( new_row != -1 );
      assert( rowidx[new_row] < 29 );
      assert(it.row() != it.col() || new_row == k);
      newA.appendRowElement(new_row, k, it.value(), rowidx[new_row]);
    }
  }

  newA.sortRows();
}

// copy Eigen matrix type system to houdini matrix type system (both A and b)
// we assume a 29 non zeros per row sparsity pattern in A
template<typename T>
void
sim_stokesSolver<T>::copySystem(
    // input
    const BlockMatrixType &A,
    const BlockVectorType &b,
    // output
    MatrixType &newA, // colmajor
    VectorType &newb) const{
  assert(!BlockMatrixType::IsRowMajor);
    #ifndef NDEBUG
    for (int k = 0; k < A.outerSize(); ++k)
    {
        typename BlockMatrixType::InnerIterator it(A,k);
        assert( it );
    }
    #endif

  auto system_size = b.size();
  newb.init(0, system_size-1);
  newA.init(system_size,/* nonzeros = */29); // assumed
  std::vector<int> rowidx(system_size, 0); // per row
  for ( int k = 0; k < system_size; ++k )
  {
    newb(k) = b[k];
    // collect non-zeros for the new matrix
    for (typename BlockMatrixType::InnerIterator it(A,k); it; ++it)
    {
      auto row = it.row();
      assert( row != -1 );
      assert( rowidx[row] < 29 );
      assert(it.row() != it.col() || row == k);
      newA.appendRowElement(row, k, it.value(), rowidx[row]);
    }
  }

  newA.sortRows();
}

template<typename T>
SolverResult
sim_stokesSolver<T>::solveBlockwiseStokes(
    const SIM_RawField & surf,
    const SIM_RawField * const* sweights,
    const SIM_RawField * const* cweights,
    const SIM_RawField & viscfield,
    const SIM_RawField & densfield,
    const SIM_RawField * const* solid_vel,
    const SIM_RawField & surf_pres,
    SIM_VectorField * valid,
    SIM_VectorField & vel) const{
  auto system_size = getNumStokesVars();
  if(!system_size)
    return NOCHANGE;
  BlockMatrixType A, H;
  BlockVectorType b(system_size);
  BlockVectorType uold = buildVelocityVector(vel, solid_vel);
  BlockVectorType ust(getNumVelocityVars());
  {
    UT_PerfMonAutoSolveEvent event(&mySolver, "Build System Blockwise");
    buildSystemBlockwise(
        A, b, H, ust, uold, surf, sweights, cweights, viscfield, densfield, solid_vel, surf_pres);
    A.prune(0, 0);
    A.makeCompressed();
  }

  if ( 0 )
  {
  if ( !isMatrixPruned(A) )
  {
    for (int k = 0; k < A.outerSize(); ++k )
    {
      typename BlockMatrixType::InnerIterator it(A,k);
      if (!it)
      {
        std::cerr << " k = " << k << std::endl;
      }
    }
    UT_VoxelArrayIteratorI vit;
    std::cerr << " u indices: " << std::endl;
    vit.setConstArray(myUIndex.field());
    for ( vit.rewind(); !vit.atEnd(); vit.advance() )
    {
      int i = vit.x(), j = vit.y(), k = vit.z();
      if ( isCollision(vit.getValue()) )
        std::cerr << i << " " << j << " " << k << ";c  ";
      else if ( isInSystem(vit.getValue()) )
        std::cerr << i << " " << j << " " << k << ";  ";
    }
    std::cerr <<  std::endl;

    std::cerr << " v indices: " << std::endl;
    vit.setConstArray(myVIndex.field());
    for ( vit.rewind(); !vit.atEnd(); vit.advance() )
    {
      int i = vit.x(), j = vit.y(), k = vit.z();
      if ( isCollision(vit.getValue()) )
        std::cerr << i << " " << j << " " << k << ";c  ";
      else if ( isInSystem(vit.getValue()) )
        std::cerr << i << " " << j << " " << k << ";  ";
    }
    std::cerr <<  std::endl;

    std::cerr << " w indices: " << std::endl;
    vit.setConstArray(myWIndex.field());
    for ( vit.rewind(); !vit.atEnd(); vit.advance() )
    {
      int i = vit.x(), j = vit.y(), k = vit.z();
      if ( isCollision(vit.getValue()) )
        std::cerr << i << " " << j << " " << k << ";c  ";
      else if ( isInSystem(vit.getValue()) )
        std::cerr << i << " " << j << " " << k << ";  ";
    }
    std::cerr <<  std::endl;

    std::cerr << " p indices: " << std::endl;
    vit.setConstArray(myCentralIndex.field());
    for ( vit.rewind(); !vit.atEnd(); vit.advance() )
    {
      int i = vit.x(), j = vit.y(), k = vit.z();
      if ( isInSystem(vit.getValue()) )
      {
        std::cerr << i << " " << j << " " << k << ";  ";
        std::cerr << "f = " << cweights[0]->field()->getValue(i,j,k) << "; ";
        std::cerr << "l = " << sweights[0]->field()->getValue(i,j,k) << "; ";
        std::cerr << "uf0 = " << cweights[4]->field()->getValue(i,j,k) << "; ";
        std::cerr << "vf0 = " << cweights[5]->field()->getValue(i,j,k) << "; ";
        std::cerr << "wf0 = " << cweights[6]->field()->getValue(i,j,k) << "; ";
        std::cerr << "ul0 = " << sweights[4]->field()->getValue(i,j,k) << "; ";
        std::cerr << "vl0 = " << sweights[5]->field()->getValue(i,j,k) << "; ";
        std::cerr << "wl0 = " << sweights[6]->field()->getValue(i,j,k) << "; ";
        std::cerr << "ub0 = " << (u_oob(i+1,j,k) || u_oob(i-1,j,k)) << "; ";
        std::cerr << "vb0 = " << (v_oob(i,j+1,k) || v_oob(i,j-1,k)) << "; ";
        std::cerr << "wb0 = " << (w_oob(i,j,k+1) || v_oob(i,j,k-1)) << "; ";

        std::cerr << "uf1 = " << cweights[4]->field()->getValue(i+1,j,k) << "; ";
        std::cerr << "vf1 = " << cweights[5]->field()->getValue(i,j+1,k) << "; ";
        std::cerr << "wf1 = " << cweights[6]->field()->getValue(i,j,k+1) << "; ";
        std::cerr << "ul1 = " << sweights[4]->field()->getValue(i+1,j,k) << "; ";
        std::cerr << "vl1 = " << sweights[5]->field()->getValue(i,j+1,k) << "; ";
        std::cerr << "wl1 = " << sweights[6]->field()->getValue(i,j,k+1) << "; ";

        std::cerr << "ub1 = " << (u_oob(i,j,k) || u_oob(i+2,j,k)) << "; ";
        std::cerr << "vb1 = " << (v_oob(i,j,k) || u_oob(i,j+2,k)) << "; ";
        std::cerr << "wb1 = " << (w_oob(i,j,k) || u_oob(i,j,k+2)) << "; ";
        std::cerr << std::endl;
      }
    }
    std::cerr << std::endl;

    return INVALID;
  }
  }

    #ifdef USE_EIGEN_SOLVER_FOR_BLOCKWISE_STOKES
    BlockVectorType x(system_size);
    auto result = solveSystemEigen(A,b,x);
    if (result == SUCCESS)
        updateVelocitiesBlockwise(uold - H*x + (1.0/dx) * ust, solid_vel, valid, vel);
    #else
    // 29 is the max non zeros per row in the stokes system
    MatrixType Ah;
    VectorType bh;
    UT_ExintArray to_original;
    if ( !isMatrixPruned(A) )
    {
        std::cerr<< "WARNING: matrix has been pruned" << std::endl;
        pruneSystem(A,b,Ah,bh,to_original);
    }
    else
    {
        copySystem(A, b, Ah, bh);
        for ( int i = 0; i < system_size; ++i )
        to_original.append(i);
    }

    VectorType xsmall(0, to_original.size()-1);

    auto result = solveSystem(Ah, bh, xsmall, mySolver.getUseOpenCL());
    if (result == SUCCESS)
    {
        UT_PerfMonAutoSolveEvent event(&mySolver, "Update Velocity");

        sim_updateVelocityParms parms(sweights, solid_vel, densfield, surf_pres,
            mySolver.getMinDensity(), mySolver.getMaxDensity());

        VectorType x(0, system_size-1);
        for ( int i = 0; i < to_original.size(); ++i )
        {
        x(to_original[i]) = xsmall(i);
        }

        for ( int axis = 0; axis < 3; ++axis )
        {
        if ( valid )
            valid->getField(axis)->makeConstant(0);
        updateVelocities(x, parms, valid, vel, axis);
        }
    }
    #endif
  return result;
}







template<typename T>
SolverResult
sim_stokesSolver<T>::solveSystemEigen(
    const BlockMatrixType &A,
    const BlockVectorType &b,
    BlockVectorType &x ) const{
  T tol = mySolver.getTolerance();

  UT_PerfMonAutoSolveEvent event(&mySolver, "Solve Blockwise System");
  auto system_size = b.size();
  if ( !system_size )
    return NOCHANGE; // nothing to do

  //PCG<T> solver;
  //solver.setTolerance(tol);
  //solver.setMaxIterations(3*system_size);
  SparseLU<T> solver;
  //SimplicialLDLT<double> solver;

  solver.compute( A );
  if ( solver.info() != Eigen::Success )
  {
    std::cout << "Compute failed: " << solver.info() << "\n";
    return FAILED;
  }

  x = solver.solve( b );
  if ( solver.info() != Eigen::Success )
  {
    std::cout << "Solve failed: ";
    switch (solver.info())
    {
      case Eigen::NumericalIssue:
        std::cout << "Numerical Issue\n"; return FAILED;
      case Eigen::NoConvergence:
        std::cout << "Did Not Converge\n"; return NOCONVERGE;
      case Eigen::InvalidInput:
        std::cout << "InvalidInput\n"; return INVALID;
      default:
        std::cout << "Unknown\n"; return FAILED;
    }

    // Print condition number if solve fails to see if there is a problem with
    // conditioning
    MatrixX<T> Adense;
    Adense = MatrixX<T>(A);
    T condition_number = Adense.inverse().norm() / Adense.norm();

    std::cerr << "k(A) = " << condition_number << std::endl;

  }

  //UT_WorkBuffer extra_info;
  //extra_info.sprintf("Iterations=%d, Error=%.6f", int(solver.iterations()), solver.error());
  //event.setExtraInfo(extra_info.buffer());

  return SUCCESS;
}



template<typename T>
void
sim_stokesSolver<T>::updateVelocitiesBlockwise(
    const BlockVectorType &unew,
    const SIM_RawField * const* solid_vel,
    SIM_VectorField * valid,
    SIM_VectorField &vel) const{
  UT_PerfMonAutoSolveEvent event(&mySolver, "Update Velocity");

  if ( valid )
    for ( int axis = 0; axis < 3; ++axis )
      valid->getField(axis)->makeConstant(0);

  // Update velocities based on the pressures and stresses determined by the solver.

  const UT_VoxelArrayF &u_solid = *solid_vel[0]->field();
  const UT_VoxelArrayF &v_solid = *solid_vel[1]->field();
  const UT_VoxelArrayF &w_solid = *solid_vel[2]->field();

  UT_VoxelArrayF &u = *vel.getField(0)->fieldNC();
  UT_VoxelArrayF &v = *vel.getField(1)->fieldNC();
  UT_VoxelArrayF &w = *vel.getField(2)->fieldNC();

  UT_VoxelArrayIteratorF vit;
  vit.setArray(&u);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    int idx = myUIndex(i,j,k);
    if (!isInSystem(idx))
    {
      vit.setValue(0);
      continue;
    }
    if ( valid )
      valid->getField(0)->fieldNC()->setValue(i,j,k,1);
    vit.setValue(isCollision(idx) ? u_solid.getValue(i,j,k) : unew[u_blk_idx(i,j,k)]);
  }

  vit.setArray(&v);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    int idx = myVIndex(i,j,k);
    if (!isInSystem(idx))
    {
      vit.setValue(0);
      continue;
    }
    if ( valid )
      valid->getField(1)->fieldNC()->setValue(i,j,k,1);
    vit.setValue(isCollision(idx) ? v_solid.getValue(i,j,k) : unew[v_blk_idx(i,j,k)]);
  }

  vit.setArray(&w);
  for ( vit.rewind(); !vit.atEnd(); vit.advance() )
  {
    int i = vit.x(), j = vit.y(), k = vit.z();
    int idx = myWIndex(i,j,k);
    if (!isInSystem(idx))
    {
      vit.setValue(0);
      continue;
    }
    if ( valid )
      valid->getField(2)->fieldNC()->setValue(i,j,k,1);
    vit.setValue(isCollision(idx) ? w_solid.getValue(i,j,k) : unew[w_blk_idx(i,j,k)]);
  }
}
