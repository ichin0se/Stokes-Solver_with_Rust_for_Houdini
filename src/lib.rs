#[no_mangle]
pub extern "C" fn rust_updateVelocities(
    i_buf: *mut i32,
    j_buf: *mut i32,
    k_buf: *mut i32,
    u_buf: *mut f32,
    valid_buf: *mut f32,
    rhox_buf: *mut f32,
    u_solid_buf: *mut f32,
    surfpres_buf: *mut f32,
    c_vol_liquid_buf: *mut f32,
    u_vol_liquid_buf: *mut f32,
    ez_vol_liquid_buf: *mut f32,
    ey_vol_liquid_buf: *mut f32,
    centralindex_buf: *mut i64,
    uindex_buf: *mut i64,
    vindex_buf: *mut i64,
    windex_buf: *mut i64,
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
    check_null!(valid_buf, "valid_buf");
    check_null!(rhox_buf, "rhox_buf");
    check_null!(u_solid_buf, "u_solid_buf");
    check_null!(surfpres_buf, "surfpres_buf");
    check_null!(c_vol_liquid_buf, "c_vol_liquid_buf");
    check_null!(u_vol_liquid_buf, "u_vol_liquid_buf");
    check_null!(ez_vol_liquid_buf, "ez_vol_liquid_buf");
    check_null!(ey_vol_liquid_buf, "ey_vol_liquid_buf");
    check_null!(centralindex_buf, "centralindex_buf");
    check_null!(uindex_buf, "uindex_buf");
    check_null!(vindex_buf, "vindex_buf");
    check_null!(windex_buf, "windex_buf");
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
    let valid_flatten = unsafe { std::slice::from_raw_parts_mut(valid_buf, len_usize * 3) };
    let rhox = unsafe { std::slice::from_raw_parts_mut(rhox_buf, len_usize) };
    let u_solid = unsafe { std::slice::from_raw_parts_mut(u_solid_buf, len_usize) };
    let surfpres = unsafe { std::slice::from_raw_parts_mut(surfpres_buf, len_usize) };
    let c_vol_liquid = unsafe { std::slice::from_raw_parts_mut(c_vol_liquid_buf, len_usize) };
    let u_vol_liquid = unsafe { std::slice::from_raw_parts_mut(u_vol_liquid_buf, len_usize) };
    let ez_vol_liquid = unsafe { std::slice::from_raw_parts_mut(ez_vol_liquid_buf, len_usize) };
    let ey_vol_liquid = unsafe { std::slice::from_raw_parts_mut(ey_vol_liquid_buf, len_usize) };
    let central_index = unsafe { std::slice::from_raw_parts_mut(centralindex_buf, len_usize) };
    let u_index = unsafe { std::slice::from_raw_parts_mut(uindex_buf, len_usize) };
    let v_index = unsafe { std::slice::from_raw_parts_mut(vindex_buf, len_usize) };
    let w_index = unsafe { std::slice::from_raw_parts_mut(windex_buf, len_usize) };
    let txx = unsafe { std::slice::from_raw_parts_mut(txx_buf, len_usize) };
    let txy = unsafe { std::slice::from_raw_parts_mut(txy_buf, len_usize) };
    let txz = unsafe { std::slice::from_raw_parts_mut(txz_buf, len_usize) };
    let p = unsafe { std::slice::from_raw_parts_mut(p_buf, len_usize) };

    println!("u element size is {}", len_usize);
    for n in 0..len_usize {
        let idx: i64 = u_index[n];
        //println!("central idx is: {}", central_index[n]);
        println!("u idx is: {}", u_index[n]);
        //println!("v idx is: {}", v_index[n]);
        //println!("w idx is: {}", w_index[n]);
        if is_collision(idx, collision_index) {
            u[n] = u_solid[n];
            valid_flatten[(3 * n) + axis as usize] = 1.0;
            println!("This Field is Collision.");
        } else if is_in_system(idx) {
            valid_flatten[(3 * n) + axis as usize] = 1.0;
            let gfp: f32 = ghost_fluid_surface_tension_pressure(
                axis,
                n,
                u_vol_liquid[n],
                surfpres,
                p,
                u_index,
                v_index,
                w_index,
                i_array,
                j_array,
                k_array,
            );
            println!("gfp: {}", gfp);
            let rho = rhox[n].clamp(minrho, maxrho);
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

            u[n] += factor
                * (c_vol_liquid[id_i_minus] * (p[id_i_minus] as f32)
                    - c_vol_liquid[n] * (p[n] as f32)
                    + ((c_vol_liquid[n] * (txx[n] as f32)
                        - c_vol_liquid[id_i_minus] * (txx[id_i_minus] as f32))
                        + (ez_vol_liquid[id_j_plus] * (txy[id_j_plus] as f32)
                            - ez_vol_liquid[n] * (txy[n] as f32))
                        + (ey_vol_liquid[id_k_plus] * (txz[id_k_plus] as f32)
                            - ey_vol_liquid[n] * (txz[n] as f32))))
                - factor * gfp;
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
    vindex: &mut [i64],
    windex: &mut [i64],
    i_array: &[i32],
    j_array: &[i32],
    k_array: &[i32],
) -> f32 {
    let uidx: i64;
    let pidx0: i64;
    let pidx1: i64 = p[idx];
    let p1: f32 = surfpres[idx];
    let p0: f32;

    if axis == 0 {
        uidx = uindex[idx];
        let id_i_minus = search_linear_index_from_vector(
            idx,
            i_array,
            i_array[idx] - 1,
            j_array,
            j_array[idx],
            k_array,
            k_array[idx],
        );
        pidx0 = p[id_i_minus];
        p0 = surfpres[id_i_minus];
    } else if axis == 1 {
        uidx = vindex[idx];
        let id_j_minus = search_linear_index_from_vector(
            idx,
            i_array,
            i_array[idx],
            j_array,
            j_array[idx] - 1,
            k_array,
            k_array[idx],
        );
        pidx0 = p[id_j_minus];
        p0 = surfpres[id_j_minus];
    } else if axis == 2 {
        uidx = windex[idx];
        let id_k_minus = search_linear_index_from_vector(
            idx,
            i_array,
            i_array[idx],
            j_array,
            j_array[idx],
            k_array,
            k_array[idx] - 1,
        );
        pidx0 = p[id_k_minus];
        p0 = surfpres[id_k_minus];
    } else {
        return 1.0;
    }

    if !is_in_system(uidx) {
        return 1.0;
    }
    if pidx0 == -2 && is_in_system(pidx1) {
        return -lerp(p1, p0, u_vol_liquid);
    } else if pidx1 == -2 && is_in_system(pidx0) {
        return lerp(p0, p1, u_vol_liquid);
    }
    1.0
}

fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}
