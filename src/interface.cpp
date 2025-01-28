#include <SYS/SYS_Types.h>
#include "SIM_Stokes.hpp"
#include "UT/UT_Interrupt.h"
#include <UT/UT_SparseMatrix.h>
#include "SIM/SIM_RawField.h"
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
#include <sys/types.h>
#include <vector>
#include <typeinfo>
using namespace Stokes;

extern "C" {
    int rust_updateVelocities(
        int*    i_buf,
        int*    j_buf,
        int*    k_buf,

        float*  u_buf,
        float*  rhox_buf,
        float*  u_solid_buf,
        float*  surfpres_buf,
        float*  c_vol_liquid_buf,
        float*  u_vol_liquid_buf,
        float*  ez_vol_liquid_buf,
        float*  ey_vol_liquid_buf,

        exint*    centralindex_buf,
        exint*    uindex_buf,
        exint*    txx_buf,
        exint*    txy_buf,
        exint*    txz_buf,
        exint*    p_buf,

        float minrho,
        float maxrho,

        int axis,
        int CollisionIndex,

        float dx,
        float dt,

        int ni,
        int nj,
        int nk,

        int total_size
    );
}
bool isInSystem(exint idx) {
    return idx >= 0;
}
// out of bounds checks
bool c_oob(int i, int j, int k, int ni, int nj, int nk)  {
    return i < 0 || i > ni-1 || j < 0 || j > nj-1 || k < 0 || k > nk-1;
}
bool tyz_oob(int i, int j, int k, int ni, int nj, int nk)  {
    return i < 0 || i > ni-1 || j < 0 || j > nj || k < 0 || k > nk;
}
bool txz_oob(int i, int j, int k, int ni, int nj, int nk) {
    return i < 0 || i > ni || j < 0 || j > nj-1 || k < 0 || k > nk;
}
bool txy_oob(int i, int j, int k, int ni, int nj, int nk)  {
    return i < 0 || i > ni || j < 0 || j > nj || k < 0 || k > nk-1;
}
bool u_oob(int i, int j, int k, int ni, int nj, int nk) {
    return i < 0 || i > ni || j < 0 || j > nj-1 || k < 0 || k > nk-1;
}
bool v_oob(int i, int j, int k, int ni, int nj, int nk) {
    return i < 0 || i > ni-1 || j < 0 || j > nj || k < 0 || k > nk-1;
}
bool w_oob(int i, int j, int k, int ni, int nj, int nk){
    return i < 0 || i > ni-1 || j < 0 || j > nj-1 || k < 0 || k > nk;
}
exint p_idx(int i, int j, int k, int ni, int nj, int nk, SIM_RawIndexField Field) {
    return c_oob(i,j,k, ni, nj, nk) ? exint(-1) : Field(i,j,k);
}
exint txx_idx(int i, int j, int k, int ni, int nj, int nk, SIM_RawIndexField Field, int NumPressure)  {
    return c_oob(i,j,k, ni, nj, nk) ? exint(-1) : (Field(i,j,k) + (isInSystem(Field(i,j,k)) ? NumPressure : 0));
}
exint tyy_idx(int i, int j, int k, int ni, int nj, int nk, SIM_RawIndexField Field, int NumPressure) {
    return c_oob(i,j,k, ni, nj, nk) ? exint(-1) : (Field(i,j,k) + (isInSystem(Field(i,j,k)) ? 2 * NumPressure : 0));
}
exint tyz_idx(int i, int j, int k, int ni, int nj, int nk, SIM_RawIndexField Field) {
    return tyz_oob(i,j,k, ni, nj, nk) ? exint(-1) : Field(i,j,k);
}
exint txz_idx(int i, int j, int k, int ni, int nj, int nk, SIM_RawIndexField Field) {
    return txz_oob(i,j,k, ni, nj, nk) ? exint(-1) : Field(i,j,k);
}
exint txy_idx(int i, int j, int k, int ni, int nj, int nk, SIM_RawIndexField Field)  {
    return txy_oob(i,j,k, ni, nj, nk) ? exint(-1) : Field(i,j,k);
}

extern "C" {
    int cxx_rust_interface( const UT_VectorT<fpreal32>& x,
                            UT_VoxelArrayIteratorF& vit,
                            SIM_VectorField* valid,
                            const UT_VoxelArrayF& density,
                            const UT_VoxelArrayF& u_solid,
                            const UT_VoxelArrayF& surfpres,
                            const UT_VoxelArrayF& c_vol_liquid,
                            const  UT_VoxelArrayF& u_vol_liquid,
                            const UT_VoxelArrayF& ez_vol_liquid,
                            const UT_VoxelArrayF& ey_vol_liquid,
                            const int& axis,
                            UT_Interrupt* boss,
                            const UT_JobInfo& info,
                            const int& ni,
                            const int& nj,
                            const int& nk,
                            const float& dt,
                            const float& dx,
                            const float& minrho,
                            const float& maxrho,
                            const int& NumPressure,
                            const SIM_RawIndexField& myUIndex,
                            const SIM_RawIndexField& myCentralIndex,
                            const SIM_RawIndexField& myTxyIndex,
                            const SIM_RawIndexField& myTxzIndex,
                            const SIM_RawIndexField& myTyzIndex,
                            const int myCollisionIndex
                            ) {

                                std::cerr << " Minrho: " << minrho << std::endl;
                                std::cerr << " Maxrho: " << maxrho << std::endl;
                                std::cerr << " Axis: " << axis << std::endl;
                                std::cerr << " NI: " << ni << std::endl;
                                std::cerr << " NJ: " << nj << std::endl;
                                std::cerr << " NK: " << nk << std::endl;
                                std::cerr << " DT: " << dt << std::endl;
                                std::cerr << " DX: " << dx << std::endl;
                                std::cerr << " MyNumPressureVars: " << NumPressure << std::endl;
                                std::cerr << " MyCollisionIndex: " << myCollisionIndex << std::endl;
                                for ( vit.rewind(); !vit.atEnd(); vit.advance() ) {
                                    int i = vit.x(), j = vit.y(), k = vit.z();
                                    std::cerr << " Density: " << density.getValue(i, j, k) << std::endl;
                                    std::cerr << " U_solid: " << u_solid.getValue(i, j, k) << std::endl;
                                    std::cerr << " Surfpres: " << surfpres.getValue(i, j, k) << std::endl;
                                    std::cerr << " C_vol_liquid: " << c_vol_liquid.getValue(i, j, k) << std::endl;
                                    std::cerr << " U_vol_liquid: " << u_vol_liquid.getValue(i, j, k) << std::endl;
                                    std::cerr << " Ez_vol_liquid: " << ez_vol_liquid.getValue(i, j, k) << std::endl;
                                    std::cerr << "Ey_vol_liquid: " << ey_vol_liquid.getValue(i, j, k) << std::endl;

                                    std::cerr << "MyUIndex: " << myUIndex(i, j, k) << std::endl;
                                    std::cerr << "MyCentralIndex: " << myCentralIndex(i, j, k) << std::endl;
                                    std::cerr << "MyTxyIndex: " << myTxyIndex(i, j, k) << std::endl;
                                    std::cerr << "MyTxzIndex: " << myTxzIndex(i, j, k) << std::endl;
                                    std::cerr << "MyTyzIndex: " << myTyzIndex(i, j, k) << std::endl;


                                }
        int len = 0;
        int flg = 0;
        for (vit.rewind(); !vit.atEnd(); vit.advance()){
            len += 1;
        }
        if(len <= 0) {
            return 0;
        }
        std::vector<float> u_buf(len+3);
        std::vector<float> valid_buf(len+3);
        std::vector<float> rhox_buf(len+3);
        std::vector<float> u_solid_buf(len+3);
        std::vector<float> surfpres_buf(len+3);
        std::vector<float> c_vol_liquid_buf(len+3);
        std::vector<float> u_vol_liquid_buf(len+3);
        std::vector<float> ez_vol_liquid_buf(len+3);
        std::vector<float> ey_vol_liquid_buf(len+3);

        std::vector<exint>   CentralIndex_buf(len+3);
        std::vector<exint>   UIndex_buf(len+3);
        std::vector<exint>   txx_buf(len+3);
        std::vector<exint>   txy_buf(len+3);
        std::vector<exint>   txz_buf(len+3);
        std::vector<exint>   p_buf(len+3);
        std::vector<int>   i_buf(len+3);
        std::vector<int>   j_buf(len+3);
        std::vector<int>   k_buf(len+3);

        UT_VoxelProbeAverage<float,-1,0,0> rhox;
        rhox.setArray(&density);

        auto txy = [&](int i, int j, int k) {
            return isInSystem(txy_idx(i, j, k, ni, nj, nk, myTxyIndex)) ? x(txy_idx(i, j, k, ni, nj, nk, myTxyIndex)) : 0;
        };
        auto txz = [&](int i, int j, int k) {
            return isInSystem(txz_idx(i, j, k, ni, nj, nk, myTxzIndex)) ? x(txz_idx(i, j, k, ni, nj, nk, myTxzIndex)) : 0;
        };
        auto txx = [&](int i, int j, int k) {
            return isInSystem(txx_idx(i, j, k, ni, nj, nk, myCentralIndex, NumPressure)) ? x(txx_idx(i, j, k, ni, nj, nk, myCentralIndex, NumPressure)) : 0;
        };
        auto p   = [&](int i, int j, int k) {
            return isInSystem(p_idx(i, j, k, ni, nj, nk, myCentralIndex))    ? x(p_idx(i, j, k, ni, nj, nk, myCentralIndex))   : 0;
        };


        if (flg == 0) {
            int lin_idx = 1;
            int i = 0;
            int j = 0;
            int k = 0;
            for ( vit.rewind(); !vit.atEnd(); vit.advance() ) {
                i = vit.x(), j = vit.y(), k = vit.z();

                    if(lin_idx == 1){
                        rhox_buf[0] = 0;
                        u_solid_buf[0] = 0;
                        surfpres_buf[0] = 0;
                        c_vol_liquid_buf[0] = 0;
                        u_vol_liquid_buf[0] = 0;
                        ez_vol_liquid_buf[0] = 0;
                        ey_vol_liquid_buf[0] = 0;
                        CentralIndex_buf[0] = 0;
                        UIndex_buf[0] = 0;
                        txy_buf[0] = 0;
                        txz_buf[0] = 0;

                        txx_buf[0] = txx(i-1, j, k);
                        p_buf[0] = p(i-1, j, k);
                        i_buf[0] = i - 1;
                        j_buf[0] = j;
                        k_buf[0] = k;
                    }

                    u_buf[lin_idx] = vit.getValue();

                    rhox.setIndex(vit);
                    rhox_buf[lin_idx] = rhox.getValue();

                    u_solid_buf[lin_idx] = u_solid.getValue(i,j,k);
                    surfpres_buf[lin_idx] = surfpres.getValue(i,j,k);

                    c_vol_liquid_buf[lin_idx] = c_vol_liquid.getValue(i,j,k);
                    u_vol_liquid_buf[lin_idx] = u_vol_liquid.getValue(i,j,k);
                    ez_vol_liquid_buf[lin_idx] = ez_vol_liquid.getValue(i,j,k);
                    ey_vol_liquid_buf[lin_idx] = ey_vol_liquid.getValue(i,j,k);

                    CentralIndex_buf[lin_idx] = (exint)myCentralIndex(i, j, k);
                    //CentralIndex_buf[lin_idx] = 314;
                    std::cerr << " CentralIndex buffer: " << CentralIndex_buf[lin_idx] << std::endl;
                    UIndex_buf[lin_idx] = (exint)myUIndex(i, j, k);
                    //UIndex_buf[lin_idx] = 314;
                    std::cerr << " UIndex buffer: " << UIndex_buf[lin_idx] << std::endl;

                    txx_buf[lin_idx] = txx(i, j, k);
                    txy_buf[lin_idx] = txy(i, j, k);
                    txz_buf[lin_idx] = txz(i, j, k);

                    p_buf[lin_idx] = p(i, j, k);

                    i_buf[lin_idx] = i;
                    j_buf[lin_idx] = j;
                    k_buf[lin_idx] = k;

                    lin_idx += 1;
                    //}
            }


            txy_buf[lin_idx] = txy(i, j + 1, k);
            txz_buf[lin_idx] = 0;
            i_buf[lin_idx] = i;
            j_buf[lin_idx] = j + 1;
            k_buf[lin_idx] = k;

            txz_buf[lin_idx + 1] = txz(i, j, k + 1);
            i_buf[lin_idx + 1] = i;
            j_buf[lin_idx + 1] = j;
            k_buf[lin_idx + 1] = k + 1;

            int ret = rust_updateVelocities(
                i_buf.data(),
                j_buf.data(),
                k_buf.data(),
                u_buf.data(),
                rhox_buf.data(),
                u_solid_buf.data(),
                surfpres_buf.data(),
                c_vol_liquid_buf.data(),
                u_vol_liquid_buf.data(),
                ez_vol_liquid_buf.data(),
                ey_vol_liquid_buf.data(),
                CentralIndex_buf.data(),
                UIndex_buf.data(),
                txx_buf.data(),
                txy_buf.data(),
                txz_buf.data(),
                p_buf.data(),
                minrho,
                maxrho,
                axis,
                myCollisionIndex,
                dx,
                dt,
                ni,
                nj,
                nk,
                (int)u_buf.size()
            );
            if (ret == 0) {
                return 0;
            }
            lin_idx = 1;
            for (vit.rewind(); !vit.atEnd(); vit.advance()){
                if (boss->opInterrupt()){
                    flg = 1;
                    break;
                } else {
                    //std::cerr << " Update Velocities by own interface." << std::endl;
                    vit.setValue(u_buf[lin_idx]);
                    valid->getField(0)->fieldNC()->setValue(vit.x(), vit.y(), vit.z(), 1);
                    lin_idx += 1;
                }
            }
        }

        if(flg == 0) {
            //std::cerr << " Done! " << std::endl;
            return 1;
        }else{
            //std::cerr << " About! " << std::endl;
            return 0;
        }
    }
}
