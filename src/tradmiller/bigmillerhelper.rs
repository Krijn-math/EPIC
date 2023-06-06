//all helper functions in the miller loop

use num_bigint::BigUint;

use crate::{bigint::*, bigconstants::BIGONE};

pub struct Point {
    pub x2: BigUint,
    pub xz: BigUint,
    pub z2: BigUint,
    pub yz: BigUint,
}

pub struct AffPoint {
    pub x: BigUint,
    pub y: BigUint,
}

pub struct Result {
    pub re: BigUint,
    pub im: BigUint,
}

pub fn point_cast(px: &BigUint, py: &BigUint) -> Point {
    let xz = (*px).clone();
    let yz = (*py).clone();
    let z2 = (*BIGONE).clone();
    let x2 = fp_sq(&xz);
    Point { x2, xz, z2, yz }
}

pub fn neg_pro(p: &Point) -> Point{
    Point { x2: p.x2.clone(), xz: p.xz.clone(), z2: p.z2.clone(), yz: fp_neg(&p.yz)}
}

pub fn neg(p: &AffPoint) -> AffPoint{
    AffPoint { x: p.x.clone(), y: fp_neg(&p.y)}
}

pub fn dbl_and_line(
    t: Point,
    mont_A: &BigUint,
) -> (
    Point,
    (BigUint, BigUint, BigUint),
) {
    let mut x2: BigUint = (t.x2).clone();
    let mut xz: BigUint = (t.xz).clone();
    let mut z2: BigUint = (t.z2).clone();
    let mut yz: BigUint = t.yz;

    let mut lx: BigUint;
    let mut ly: BigUint;
    let mut l0: BigUint;

    let mut xx2: BigUint;

    xx2 = fp_add(&yz, &yz);
    ly = fp_sq(&xx2);
    l0 = fp_sub(&x2, &z2);
    let v0: BigUint = fp_sq(&l0);
    l0 = fp_mul(&xx2, &l0);
    lx = fp_mul(&xz, &l0);
    xx2 = fp_mul(&yz, &ly);
    lx = fp_add(&xx2, &lx);
    yz = fp_add(&x2, &z2);
    yz = fp_mul(mont_A, &yz);
    xx2 = fp_add(&xz, &xz);
    yz = fp_add(&xx2, &yz);
    yz = fp_add(&xx2, &yz);
    yz = fp_mul(&xx2, &yz);

    xx2 = fp_sq(&v0);
    let t0: BigUint = fp_sq(&l0);
    z2 = fp_sq(&ly);
    yz = fp_add(&v0, &yz);
    yz = fp_mul(&l0, &yz);

    ly = fp_mul(&xz, &ly);
    l0 = fp_mul(&x2, &l0);

    x2 = xx2;
    xz = t0;

    (Point {x2, xz, z2, yz}, (lx, ly, l0))
}

pub fn give_line(
    lx: &BigUint,
    ly: &BigUint,
    l0: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
) -> (BigUint, BigUint) {
    let mut resre = fp_mul(lx, rx);
    let resim = fp_mul(ly, ry);

    resre = fp_sub(l0, &resre);

    (resre, resim)
}

pub fn square_and_absorb(
    f: Result,
    ellxre: &BigUint,
    ellxim: &BigUint,
) -> Result {
    let (mut re, mut im) = fp2sqr(&f.re, &f.im);
    (re, im) = fp2mul(&re, &im, ellxre, ellxim);
    Result{ re, im}
}

pub fn add_and_line(
    t: Point,
    p: &AffPoint,
    rx: &BigUint,
    mont_A: &BigUint,
) -> (Point, (BigUint, BigUint)) {
    let mut lambdax: BigUint;
    let mut lambdaz: BigUint;
    let mut ellxre: BigUint;

    lambdax = fp_mul(&p.y, &t.z2);
    lambdax = fp_sub(&lambdax, &t.yz);

    lambdaz = fp_mul(&p.x, &t.z2); //additional trick is that px can be low so implemented as adds
    lambdaz = fp_sub(&lambdaz, &t.xz);

    let lambdaz2 = fp_sq(&lambdaz); //tmp1 = pow(lambdaz, 2);
    let tylambdaz2 = fp_mul(&t.yz, &lambdaz2);

    //Get ell = line(t, p, r) = (ry - ty - lambda(rx - tx)) / (rx + tx + px - lambda*lambda + a)
    //note here that rx and px are related for specific points as sampled by us, namely rx = -px

    ellxre = fp_mul(rx, &t.z2);
    ellxre = fp_sub(&t.xz, &ellxre);

    ellxre = fp_mul(&ellxre, &lambdaz);
    ellxre = fp_mul(&ellxre, &lambdax);
    ellxre = fp_sub(&ellxre, &tylambdaz2);

    //it feels like we can improve this by using the rep (XX : XZ : ZZ : YZ)
    //right now we just use the (XZ : ZZ : YZ) part and fix it in the end with v
    //perhaps all of this can also be combined well with ellxre

    let tzlambdaz2 = fp_mul(&t.z2, &lambdaz2); //tmp2 = tz * tmp1;
    let lambdax2 = fp_sq(&lambdax); //tmp3 = pow(lambdax, 2);
    let txlambdaz2 = fp_mul(&t.xz, &lambdaz2); //tmp4 = tx * tmp1; //used in step 1.b and step 3
    let lambdaz3 = fp_mul(&lambdaz2, &lambdaz); //tmp5 = lambdaz * tmp1;
    let tzlambdax2 = fp_mul(&t.z2, &lambdax2); //tmp6 = tz * tmp3;

    let mut newtx = fp_add(&p.x, mont_A);
    newtx = fp_mul(&newtx, &fp_neg(&tzlambdaz2));
    newtx = fp_sub(&newtx, &txlambdaz2);
    newtx = fp_add(&newtx, &tzlambdax2);

    let tmp = fp_mul(&lambdaz3, &t.yz);
    let mut newty = fp_sub(&txlambdaz2, &newtx);
    newty = fp_mul(&newty, &lambdax);
    newty = fp_sub(&newty, &tmp);

    let x2 = fp_sq(&newtx);
    let xz = fp_mul(&tzlambdaz2, &newtx);
    let z2 = fp_sq(&tzlambdaz2);
    let mut yz = fp_mul(&lambdaz, &newty);
    yz = fp_mul(&t.z2, &yz);

    (Point {x2, xz, z2, yz}, (ellxre, tzlambdaz2))
}

pub fn give_line_bit(lambdax: &BigUint, lambday: &BigUint, ry: &BigUint) -> (BigUint, BigUint) {
    let resim = fp_mul(ry, lambday);

    ((*lambdax).clone(), resim)
}

pub fn absorb(
    f: Result,
    ellxre: &BigUint,
    ellxim: &BigUint,
) -> Result {
    let (re, im) = fp2mul(&f.re, &f.im, ellxre, ellxim);
    Result{ re, im }
}