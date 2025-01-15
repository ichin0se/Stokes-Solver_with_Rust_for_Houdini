// src/lib.rs

use std::os::raw::{c_double, c_float, c_int};

//
// 1) Realトレイト & impl
//
pub trait Real: Copy + Default + 'static {
    fn from_f32(x: f32) -> Self;
    fn from_f64(x: f64) -> Self;
    fn add_one(x: Self) -> Self;
    fn add_ten(x: Self) -> Self;
}

/// f32 用実装
impl Real for f32 {
    fn from_f32(x: f32) -> Self {
        x
    }
    fn from_f64(x: f64) -> Self {
        x as f32
    }
    fn add_one(x: Self) -> Self {
        x + 1.0
    }
    fn add_ten(x: Self) -> Self {
        x + 10.0
    }
}

/// f64 用実装
impl Real for f64 {
    fn from_f32(x: f32) -> Self {
        x as f64
    }
    fn from_f64(x: f64) -> Self {
        x
    }
    fn add_one(x: Self) -> Self {
        x + 1.0
    }
    fn add_ten(x: Self) -> Self {
        x + 10.0
    }
}

//
// 2) SolverResult
//
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverResult {
    Success = 1,
    Fail = 0,
}

//
// 3) SimStokesSolver<T: Real> 構造体
//
pub struct SimStokesSolver<T: Real> {
    pub nx: c_int,
    pub ny: c_int,
    pub nz: c_int,
    pub dx: T,
    pub dt: T,
}

impl<T: Real> SimStokesSolver<T> {
    pub fn new(nx: c_int, ny: c_int, nz: c_int, dx: T, dt: T) -> Self {
        SimStokesSolver { nx, ny, nz, dx, dt }
    }

    pub fn classify_and_build_indices(&mut self) {
        eprintln!(
            "(Rust) classify_and_build_indices: (nx,ny,nz)=({},{},{})",
            self.nx, self.ny, self.nz
        );
        // ここではダミー
    }

    /// solve: velocity, surface, collision の配列を受け取り
    /// - velocity を +1 or +10 (f32/f64区別)
    /// - surface/collision が空なら失敗 というダミー例
    pub fn solve(&mut self, velocity: &mut [T], surface: &[T], collision: &[T]) -> SolverResult {
        if velocity.is_empty() || surface.is_empty() {
            return SolverResult::Fail;
        }

        // f32 → +1, f64 → +10 など区別してみる（単なる例）
        for v in velocity {
            // T::add_one / T::add_ten
            *v = if std::any::type_name::<T>() == "f32" {
                T::add_one(*v)
            } else {
                T::add_ten(*v)
            };
        }

        if collision.is_empty() {
            return SolverResult::Fail;
        }

        SolverResult::Success
    }
}

//
// 4) extern "C" で f32/f64 専用関数をエクスポート
//    -> C++ 側から呼び出せるようにする
//

/// f32 create
#[no_mangle]
pub extern "C" fn stokes_solver_create_f32(
    nx: c_int,
    ny: c_int,
    nz: c_int,
    dx: c_float,
    dt: c_float,
) -> *mut SimStokesSolver<f32> {
    let solver = SimStokesSolver::new(nx, ny, nz, f32::from_f32(dx), f32::from_f32(dt));
    Box::into_raw(Box::new(solver))
}

/// f32 destroy
#[no_mangle]
pub extern "C" fn stokes_solver_destroy_f32(ptr: *mut SimStokesSolver<f32>) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Box::from_raw(ptr); // drop
    }
}

/// f32 classify
#[no_mangle]
pub extern "C" fn stokes_solver_classify_f32(ptr: *mut SimStokesSolver<f32>) -> c_int {
    if ptr.is_null() {
        return SolverResult::Fail as c_int;
    }
    let solver = unsafe { &mut *ptr };
    solver.classify_and_build_indices();
    SolverResult::Success as c_int
}

/// f32 solve
#[no_mangle]
pub extern "C" fn stokes_solver_solve_f32(
    solver_ptr: *mut SimStokesSolver<f32>,
    velocity_buf: *mut c_float,
    velocity_len: c_int,
    surface_buf: *const c_float,
    surface_len: c_int,
    collision_buf: *const c_float,
    collision_len: c_int,
) -> c_int {
    if solver_ptr.is_null() {
        return SolverResult::Fail as c_int;
    }
    let solver = unsafe { &mut *solver_ptr };

    // スライスへ変換
    let vel_slice = unsafe { std::slice::from_raw_parts_mut(velocity_buf, velocity_len as usize) };
    let surf_slice = unsafe { std::slice::from_raw_parts(surface_buf, surface_len as usize) };
    let coll_slice = unsafe { std::slice::from_raw_parts(collision_buf, collision_len as usize) };

    let res = solver.solve(vel_slice, surf_slice, coll_slice);
    res as c_int
}

// ---- f64 同様に ----

#[no_mangle]
pub extern "C" fn stokes_solver_create_f64(
    nx: c_int,
    ny: c_int,
    nz: c_int,
    dx: c_double,
    dt: c_double,
) -> *mut SimStokesSolver<f64> {
    let solver = SimStokesSolver::new(nx, ny, nz, f64::from_f64(dx), f64::from_f64(dt));
    Box::into_raw(Box::new(solver))
}

#[no_mangle]
pub extern "C" fn stokes_solver_destroy_f64(ptr: *mut SimStokesSolver<f64>) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub extern "C" fn stokes_solver_classify_f64(ptr: *mut SimStokesSolver<f64>) -> c_int {
    if ptr.is_null() {
        return SolverResult::Fail as c_int;
    }
    let solver = unsafe { &mut *ptr };
    solver.classify_and_build_indices();
    SolverResult::Success as c_int
}

#[no_mangle]
pub extern "C" fn stokes_solver_solve_f64(
    solver_ptr: *mut SimStokesSolver<f64>,
    velocity_buf: *mut c_double,
    velocity_len: c_int,
    surface_buf: *const c_double,
    surface_len: c_int,
    collision_buf: *const c_double,
    collision_len: c_int,
) -> c_int {
    if solver_ptr.is_null() {
        return SolverResult::Fail as c_int;
    }
    let solver = unsafe { &mut *solver_ptr };

    let vel_slice = unsafe { std::slice::from_raw_parts_mut(velocity_buf, velocity_len as usize) };
    let surf_slice = unsafe { std::slice::from_raw_parts(surface_buf, surface_len as usize) };
    let coll_slice = unsafe { std::slice::from_raw_parts(collision_buf, collision_len as usize) };

    let res = solver.solve(vel_slice, surf_slice, coll_slice);
    res as c_int
}

//
// 5) テスト (複数 mod でOK)
//
#[cfg(test)]
mod tests {
    use super::*;
    use std::mem;

    #[test]
    fn test_f32_success() {
        let ptr = stokes_solver_create_f32(10, 20, 30, 1.0, 0.1);
        assert!(!ptr.is_null());
        let ret_cls = stokes_solver_classify_f32(ptr);
        assert_eq!(ret_cls, SolverResult::Success as i32);

        let mut velocity = vec![0f32; 5];
        let surface = vec![1f32; 5];
        let collision = vec![2f32; 5];

        let ret_solve = stokes_solver_solve_f32(
            ptr,
            velocity.as_mut_ptr(),
            velocity.len() as i32,
            surface.as_ptr(),
            surface.len() as i32,
            collision.as_ptr(),
            collision.len() as i32,
        );
        assert_eq!(ret_solve, SolverResult::Success as i32);

        // f32 => +1
        for &val in &velocity {
            assert_eq!(val, 1.0);
        }

        stokes_solver_destroy_f32(ptr);
    }

    #[test]
    fn test_f32_fail() {
        let ptr = stokes_solver_create_f32(10, 20, 30, 1.0, 0.1);
        let mut velocity = vec![0f32; 5];
        let surface = vec![1f32; 5];
        let collision = Vec::<f32>::new(); // empty => fail

        let ret_solve = stokes_solver_solve_f32(
            ptr,
            velocity.as_mut_ptr(),
            velocity.len() as i32,
            surface.as_ptr(),
            surface.len() as i32,
            collision.as_ptr(),
            collision.len() as i32,
        );
        assert_eq!(ret_solve, SolverResult::Fail as i32);

        stokes_solver_destroy_f32(ptr);
    }

    #[test]
    fn test_f64_success() {
        let ptr = stokes_solver_create_f64(10, 20, 30, 1.0, 0.1);
        assert!(!ptr.is_null());
        let ret_cls = stokes_solver_classify_f64(ptr);
        assert_eq!(ret_cls, 1);

        let mut velocity = vec![0f64; 5];
        let surface = vec![1f64; 5];
        let collision = vec![2f64; 5];

        let ret_solve = stokes_solver_solve_f64(
            ptr,
            velocity.as_mut_ptr(),
            velocity.len() as i32,
            surface.as_ptr(),
            surface.len() as i32,
            collision.as_ptr(),
            collision.len() as i32,
        );
        assert_eq!(ret_solve, 1);

        // f64 => +10
        for &val in &velocity {
            assert_eq!(val, 10.0);
        }

        stokes_solver_destroy_f64(ptr);
    }

    #[test]
    fn test_f64_fail() {
        let ptr = stokes_solver_create_f64(10, 20, 30, 1.0, 0.1);
        let mut velocity = vec![0f64; 5];
        let surface = vec![1f64; 5];
        let collision = Vec::<f64>::new();

        let ret_solve = stokes_solver_solve_f64(
            ptr,
            velocity.as_mut_ptr(),
            velocity.len() as i32,
            surface.as_ptr(),
            surface.len() as i32,
            collision.as_ptr(),
            collision.len() as i32,
        );
        assert_eq!(ret_solve, 0);

        stokes_solver_destroy_f64(ptr);
    }
}
