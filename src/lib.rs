#[no_mangle]
pub unsafe extern "C" fn rust_updateVelocities(
    i_buf: *mut i32,
    j_buf: *mut i32,
    k_buf: *mut i32,

    u_buf: *mut f32,
    rhox_buf: *mut f32,
    u_solid_buf: *mut f32,
    surfpres_buf: *mut f32,
    c_vol_liquid_buf: *mut f32,
    u_vol_liquid_buf: *mut f32,
    ez_vol_liquid_buf: *mut f32,
    ey_vol_liquid_buf: *mut f32,

    centralindex_buf: *mut i64,
    uindex_buf: *mut i64,
    txx_buf: *mut i64,
    txy_buf: *mut i64,
    txz_buf: *mut i64,
    p_buf: *mut i64,

    minrho: f32,
    maxrho: f32,

    axis: i32,
    collision_index: i32,

    dx: f32,
    dt: f32,

    ni: i32,
    nj: i32,
    nk: i32,

    total_size: i32,
) -> i32 {
    println!("Processing Rust Function!");
    println!("total size: {}", total_size);
    if total_size <= 0 {
        println!("Buffer Size if 0. Rust Functon Done!");
        return 0;
    }

    let mut is_null = false;
    macro_rules! check_null {
        ($ptr:expr, $name:expr) => {
            if $ptr.is_null() {
                println!("Pointer `{}` is NULL!", $name);
                is_null = true;
            }
        };
    }

    check_null!(i_buf, "i_buf");
    check_null!(j_buf, "j_buf");
    check_null!(k_buf, "k_buf");
    check_null!(u_buf, "u_buf");
    check_null!(rhox_buf, "rhox_buf");
    check_null!(u_solid_buf, "u_solid_buf");
    check_null!(surfpres_buf, "surfpres_buf");
    check_null!(c_vol_liquid_buf, "c_vol_liquid_buf");
    check_null!(u_vol_liquid_buf, "u_vol_liquid_buf");
    check_null!(ez_vol_liquid_buf, "ez_vol_liquid_buf");
    check_null!(ey_vol_liquid_buf, "ey_vol_liquid_buf");
    check_null!(centralindex_buf, "centralindex_buf");
    check_null!(uindex_buf, "uindex_buf");
    check_null!(txx_buf, "txx_buf");
    check_null!(txy_buf, "txy_buf");
    check_null!(txz_buf, "txz_buf");
    check_null!(p_buf, "p_buf");

    if is_null {
        println!("At least one pointer was NULL. Rust function done!");
        return 0;
    }

    let len_usize = total_size as usize;
    let i_array = unsafe { std::slice::from_raw_parts_mut(i_buf, len_usize) };
    let j_array = unsafe { std::slice::from_raw_parts_mut(j_buf, len_usize) };
    let k_array = unsafe { std::slice::from_raw_parts_mut(k_buf, len_usize) };
    let u = unsafe { std::slice::from_raw_parts_mut(u_buf, len_usize) };

    let rhox = unsafe { std::slice::from_raw_parts_mut(rhox_buf, len_usize) };
    let u_solid = unsafe { std::slice::from_raw_parts_mut(u_solid_buf, len_usize) };
    let surfpres = unsafe { std::slice::from_raw_parts_mut(surfpres_buf, len_usize) };
    let c_vol_liquid = unsafe { std::slice::from_raw_parts_mut(c_vol_liquid_buf, len_usize) };
    let u_vol_liquid = unsafe { std::slice::from_raw_parts_mut(u_vol_liquid_buf, len_usize) };
    let ez_vol_liquid = unsafe { std::slice::from_raw_parts_mut(ez_vol_liquid_buf, len_usize) };
    let ey_vol_liquid = unsafe { std::slice::from_raw_parts_mut(ey_vol_liquid_buf, len_usize) };
    let central_index = unsafe { std::slice::from_raw_parts_mut(centralindex_buf, len_usize) };
    let u_index = unsafe { std::slice::from_raw_parts_mut(uindex_buf, len_usize) };
    let txx = unsafe { std::slice::from_raw_parts_mut(txx_buf, len_usize) };
    let txy = unsafe { std::slice::from_raw_parts_mut(txy_buf, len_usize) };
    let txz = unsafe { std::slice::from_raw_parts_mut(txz_buf, len_usize) };
    let p = unsafe { std::slice::from_raw_parts_mut(p_buf, len_usize) };

    for n in 1..len_usize - 3 {
        let idx: i64 = u_index[n];
        if is_collision(idx, collision_index) {
            u[n] = u_solid[n];
        } else if is_in_system(idx) {
            let gfp: f32 = ghost_fluid_surface_tension_pressure(
                axis,
                n,
                u_vol_liquid[n],
                surfpres,
                p,
                u_index,
                central_index,
                i_array,
                j_array,
                k_array,
                ni,
                nj,
                nk,
            );

            let rho = rhox[n].clamp(minrho, maxrho);
            println!("uindex: {}", u_index[n]);
            println!("central_index: {}", central_index[n]);
            println!("u_vol_liquid: {}", u_vol_liquid[n]);
            let factor: f32 = dt / (dx * rho * u_vol_liquid[n]);

            // The meaning of id_i_minus is a linear index `n` representing the index `(i-1, j, k)`
            let id_i_minus = search_linear_index_from_vector(
                n,
                i_array,
                i_array[n] - 1,
                j_array,
                j_array[n],
                k_array,
                k_array[n],
            );

            // The meaning of id_j_plus is a linear index `n` representing the index `(i, j+1, k)`
            let id_j_plus = search_linear_index_from_vector(
                n,
                i_array,
                i_array[n],
                j_array,
                j_array[n] + 1,
                k_array,
                k_array[n],
            );

            // The meaning of id_k_plus is a linear index `n` representing the index `(i, j, k+1)`
            let id_k_plus = search_linear_index_from_vector(
                n,
                i_array,
                i_array[n],
                j_array,
                j_array[n],
                k_array,
                k_array[n] + 1,
            );

            u[n] = (u[n]
                + factor
                    * (c_vol_liquid[id_i_minus] * (p[id_i_minus] as f32)
                        - c_vol_liquid[n] * (p[n] as f32)
                        + ((c_vol_liquid[n] * (txx[n] as f32)
                            - c_vol_liquid[id_i_minus] * (txx[id_i_minus] as f32))
                            + (ez_vol_liquid[id_j_plus] * (txy[id_j_plus] as f32)
                                - ez_vol_liquid[n] * (txy[n] as f32))
                            + (ey_vol_liquid[id_k_plus] * (txz[id_k_plus] as f32)
                                - ey_vol_liquid[n] * (txz[n] as f32))))
                - factor * gfp);
            // println!("u: {}", u[n]);
            //println!("Compute Correct Velocity Field!");
        } else {
            //println!("Fill Velocity Field by Zero.");
            u[n] = 0 as f32;
        }
    }
    1
}

/// Search for linear index `n` mapping to index `(i, j, k)`
fn search_linear_index_from_vector(
    default_idx: usize,
    i_array: &[i32],
    i: i32,
    j_array: &[i32],
    j: i32,
    k_array: &[i32],
    k: i32,
) -> usize {
    let mut index: i32 = -1;
    for n in 0..i_array.len() {
        if i_array[n] == i && j_array[n] == j && k_array[n] == k {
            index = n as i32;
            break;
        }
    }
    if index == -1 {
        index = default_idx as i32;
    }
    index as usize
}
fn is_in_system(idx: i64) -> bool {
    idx >= 0
}

fn is_collision(idx: i64, collision_index: i32) -> bool {
    idx == -3 || idx >= collision_index as i64
}

fn ghost_fluid_surface_tension_pressure(
    axis: i32,
    idx: usize,
    u_vol_liquid: f32,
    surfpres: &mut [f32],
    p: &mut [i64],
    uindex: &mut [i64],
    centralindex: &mut [i64],
    i_array: &[i32],
    j_array: &[i32],
    k_array: &[i32],
    ni: i32,
    nj: i32,
    nk: i32,
) -> f32 {
    let uidx: i64;
    let pidx0: i64;
    let pidx1: i64 = centralindex[idx];
    let p1: f32 = surfpres[idx];
    let p0: f32;

    uidx = uindex[idx];
    pidx0 = p_idx(
        i_array[idx] - 1,
        j_array[idx],
        k_array[idx],
        ni,
        nj,
        nk,
        idx,
        i_array,
        j_array,
        k_array,
        p,
    );

    let (ic, jc, kc) = clamp_ijk_to_array(
        i_array[idx] - 1,
        j_array[idx],
        k_array[idx],
        i_array,
        j_array,
        k_array,
    );
    let id_i_minus = search_linear_index_from_vector(idx, i_array, ic, j_array, jc, k_array, kc);
    p0 = surfpres[id_i_minus];

    //println!("uidx: {}", uidx);
    //println!("pidx0 : {}", pidx0);
    //println!("pidx1: {}", pidx1);
    //println!("p0: {}", p0);
    //println!("p1: {}", p1);
    // println!("u_vol_liquid: {}", u_vol_liquid);
    if !is_in_system(uidx) {
        return 0.0;
    }
    if pidx0 == -2 && is_in_system(pidx1) {
        return -lerp(p1, p0, u_vol_liquid);
    } else if pidx1 == -2 && is_in_system(pidx0) {
        return lerp(p0, p1, u_vol_liquid);
    }
    0.0
}

fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}
#[inline]
pub fn c_oob(i: i32, j: i32, k: i32, ni: i32, nj: i32, nk: i32) -> bool {
    i < 0 || i > (ni - 1) || j < 0 || j > (nj - 1) || k < 0 || k > (nk - 1)
}
pub fn p_idx(
    i: i32,
    j: i32,
    k: i32,
    ni: i32,
    nj: i32,
    nk: i32,
    default_idx: usize,
    i_array: &[i32],
    j_array: &[i32],
    k_array: &[i32],
    data: &mut [i64],
) -> i64 {
    if c_oob(i, j, k, ni, nj, nk) {
        -1.0 as i64
    } else {
        let index =
            search_linear_index_from_vector(default_idx, i_array, i, j_array, j, k_array, k);
        data[index]
    }
}
fn clamp_ijk_to_array(
    i: i32,
    j: i32,
    k: i32,
    i_array: &[i32],
    j_array: &[i32],
    k_array: &[i32],
) -> (i32, i32, i32) {
    if i_array.is_empty() || j_array.is_empty() || k_array.is_empty() {
        return (i, j, k);
    }

    let i_min = *i_array.iter().min().unwrap();
    let i_max = *i_array.iter().max().unwrap();
    let j_min = *j_array.iter().min().unwrap();
    let j_max = *j_array.iter().max().unwrap();
    let k_min = *k_array.iter().min().unwrap();
    let k_max = *k_array.iter().max().unwrap();

    let i_clamped = i.clamp(i_min, i_max);
    let j_clamped = j.clamp(j_min, j_max);
    let k_clamped = k.clamp(k_min, k_max);

    (i_clamped, j_clamped, k_clamped)
}
