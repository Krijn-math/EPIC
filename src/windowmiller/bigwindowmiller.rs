//here, I've tried to implement the approach from 2016/963 to optimize the arith for T <- T + T and (fxre, fxim)

use crate::bigconstants::*;
use crate::bigint::*;
use crate::tradmiller::bigmillerhelper::*;
use crate::windowmiller::bigger_window::*;
use num_bigint::BigUint;

pub fn double(
    t: Point,
    f: Result,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (Point, Result) {
    //check if t is inf or two
    if is_inf(&t.xz, &t.yz, &t.z2) | is_two(&t.xz, &t.yz, &t.z2) {
        //square f, set t to INF
        let (re, im) = fp2sqr(&f.re, &f.im);
        let x2 = t.x2.clone();
        let (xz, yz, z2) = ((*BIGZERO).clone(), (*BIGONE).clone(), (*BIGZERO).clone());
        let t_res: Point = Point { x2, xz, yz, z2 };
        let f_res: Result = Result { re, im };
        return (t_res, f_res);
    } else {
        //double T, get line coeffs
        let (t_res, (lambdax, lambday, lambdaz)) = dbl_and_line(t, mont_A);

        //evaluate line in rx, ry
        let (ellxre, ellxim) = give_line(&lambdax, &lambday, &lambdaz, rx, ry);

        //square f and multiply by evaluation
        let f_res = square_and_absorb(f, &ellxre, &ellxim);
        return (t_res, f_res);
    }
}

pub fn n_double(
    n: usize,
    t: Point,
    f: Result,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (Point, Result) {
    let mut t_res = t;
    let mut f_res = f;

    for _i in 0..n {
        (t_res, f_res) = double(t_res, f_res, &rx, &ry, mont_A);
    }

    (t_res, f_res)
}

pub fn add(
    t: Point,
    p: &AffPoint,
    f: Result,
    fadd: &Result,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (Point, Result) {
    if is_inf(&t.xz, &t.yz, &t.z2) {
        return (t, f);
    }

    //if t is two, we do nothing will ell or f
    if is_two(&t.xz, &t.yz, &t.z2) {
        //add P to T, get line coeffs, here tx2 is not actually computed
        let (t_res, (lambdax, lambday)) = add_and_line(t, &p, rx, mont_A);
        let re = (f.re).clone();
        let im = (f.im).clone();
        return (t_res, Result { re, im });
    }

    //add P to T, get line coeffs, here tx2 is not actually computed
    let (res_t, (lambdax, lambday)) = add_and_line(t, &p, rx, mont_A);
    let (ellxre, ellxim) = give_line_bit(&lambdax, &lambday, ry);
    let res_f = absorb(f, &ellxre, &ellxim);
    let (resre, resim) = fp2mul(&res_f.re, &res_f.im, &fadd.re, &fadd.im);

    (res_t, Result{ re: resre, im: resim})
}

pub fn subtract(
    t: Point,
    p: &AffPoint,
    f: Result,
    fsub: &Result,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (Point, Result) {
    add(t, &neg(&p), f, fsub, &rx, &ry, mont_A)
}

pub fn windowmiller(
    _n: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> Result {
    //cleanmiller but implemented using a window and NAF, using the fact that we can easily add or substract

    let t: Point = point_cast(&px, &py);
    let mut f: Result = Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() };

    let p: AffPoint = AffPoint {
        x: (*px).clone(),
        y: (*py).clone(),
    };

    //precomputing the window
    let (p1, p3, p5, p7, p9, p11, p13, p15, p17, p19, p21, f1, f1min, f3, f3min, f5, f5min, f7, f7min, f9, f9min, f11, f11min, f13, f13min, f15, f15min, f17, f17min, f19, f19min, f21, f21min) =
        precompute(&p, &rx, &ry, &mont_A);

    f = manual_window(t, f, &rx, &ry, &mont_A, &p1, &p3, &p5, &p7, &p9, &p11, &p13, &p15, &p17, &p19, &p21, &f1, &f1min, &f3, &f3min, &f5, &f5min, &f7, &f7min, &f9, &f9min, &f11, &f11min, &f13, &f13min, &f15, &f15min, &f17, &f17min, &f19, &f19min, &f21, &f21min);
    let (re, im) = fp2sqr(&f.re, &f.im);        //last squaring instead of a doubling of T and sq f
    Result {re, im}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ezero_test() {
        //this is a verified example with magma for E0
        let mont_A: BigUint = (*BIGZERO).clone();

        let tmx: BigUint = BigUint::parse_bytes(b"7", 10).unwrap();
        let tmy: BigUint = BigUint::parse_bytes(b"2D6D8560086EF2C1C9ECB70B54BE94C5D4089051E5CD25E1C812AAFB4D0B876019A9BDC1259CCED43FC83FE5B4D277080BAB880649618832D87643A7111497D4", 16).unwrap();

        let tpx: BigUint = fp_neg(&tmx); //EXTRASPECIAL FORM
        let tpy: BigUint = tmy.clone(); //E0 form

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let f = windowmiller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&f.re, &f.im, &to_red);

        //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"43E0EDF3667B03F03E05B72AC7D87032CE217C9837B6A28E0BEB265AE92D37A7B694DEA79B3D06711631B3CF92B2E5F8FE537F630CACC01D98562E2A8DDB2AA", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"20701C645F34AE64D163CB21503A8FA0D3B4521E890EEA504ED152C53709DCCD234AFC78A7B3308E1A48E0339FC75E9151A168D12A9501B3024AE552B9430804", 16).unwrap();
        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }

    #[test]
    fn ethree_test() {
        //this is a verified example with magma for E0
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();

        let tmx: BigUint = BigUint::parse_bytes(b"2", 10).unwrap();
        let tmy: BigUint = BigUint::parse_bytes(b"15871D66F2CE52EE533800A4716AD8DEB15022BE55E48E041BB09FBD31BE148E164FFCF7B10E76E4AAC5946FD3C072765C14FEEFDE36893CB78B23161BFB6D1A", 16).unwrap();

        let tpx: BigUint = fp_neg(&tmx); //EXTRASPECIAL FORM
        let tpy: BigUint = BigUint::parse_bytes(b"12C9163DF1EC9421CF06E661549E2C9069143312489A7C7AD4F41BF2C2F8FB90954C96691F3691649FCB4FD89C6E398462F72D441778945E39D8414C654AB56D", 16).unwrap();

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let f = windowmiller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&f.re, &f.im, &to_red);

        //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"7230A53FA464DF7559C7716823F50D114A2A96069A89F0027C01837A5BF7912642D0897BE621FADC7134A98E787C135501F462C1AEFE5DD1DDA8ED8BBE42865", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"2892CA051611D66B7A0DE4470C028535C04A57870FFA9C5FAC75644774EDE7AB2C9728E62D54CC77E33D7DC61088F4D6C873D167D0651F620C3904AA374C8C69", 16).unwrap();

        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }
}
