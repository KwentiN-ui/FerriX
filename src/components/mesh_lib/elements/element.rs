/// <https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node194.html>
#[derive(Debug)]
pub enum Element {
    // General 3D-Solids
    /// 4-node linear tetrahedral element
    C3D4(usize, [usize; 4]),
    /// 6-node linear triangular prism element
    C3D6(usize, [usize; 6]),
    /// 3D 20-node quadratic isoparametric element
    C3D20(usize, [usize; 20]),
}
