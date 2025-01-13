#[no_mangle]
pub extern "C" fn stokes_solver(dstbuf: *mut f32, srcbuf: *mut f32, len: i32) -> i32 {
    if dstbuf.is_null() || srcbuf.is_null() || len <= 0 {
        return 0;
    }

    1
}
