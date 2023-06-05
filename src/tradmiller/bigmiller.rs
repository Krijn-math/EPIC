//here, I've tried to implement the approach from 2016/963 to optimize the arith for T <- T + T and (fxre, fxim)


use crate::bigconstants::*;
use crate::bigint::*;
use num_bigint::BigUint;

pub fn dbl_and_line(
    x2: &BigUint,
    xz: &BigUint,
    z2: &BigUint,
    yz: &BigUint,
    mont_A: &BigUint,
) -> (
    (BigUint, BigUint, BigUint, BigUint),
    (BigUint, BigUint, BigUint),
) {
    let mut vx2: BigUint = (*x2).clone();
    let mut vxz: BigUint = (*xz).clone();
    let mut vz2: BigUint = (*z2).clone();
    let mut vyz: BigUint = (*yz).clone();

    let mut lx: BigUint;
    let mut ly: BigUint;
    let mut l0: BigUint;
    let v0: BigUint;

    let mut xx2: BigUint;
    let t0: BigUint;

    xx2 = fp_add(&vyz, &vyz);
    ly = fp_sq(&xx2);
    l0 = fp_sub(&vx2, &vz2);
    v0 = fp_sq(&l0);
    l0 = fp_mul(&xx2, &l0);
    lx = fp_mul(&vxz, &l0);
    xx2 = fp_mul(&vyz, &ly);
    lx = fp_add(&xx2, &lx);
    vyz = fp_add(&vx2, &vz2);
    vyz = fp_mul(&mont_A, &vyz);
    xx2 = fp_add(&vxz, &vxz);
    vyz = fp_add(&xx2, &vyz);
    vyz = fp_add(&xx2, &vyz);
    vyz = fp_mul(&xx2, &vyz);

    xx2 = fp_sq(&v0);
    t0 = fp_sq(&l0);
    vz2 = fp_sq(&ly);
    vyz = fp_add(&v0, &vyz);
    vyz = fp_mul(&l0, &vyz);

    ly = fp_mul(&vxz, &ly);
    l0 = fp_mul(&vx2, &l0);

    vx2 = xx2;
    vxz = t0;

    ((vx2, vxz, vz2, vyz), (lx, ly, l0))
}

pub fn give_line(
    lx: &BigUint,
    ly: &BigUint,
    l0: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
) -> (BigUint, BigUint) {
    let mut resre = fp_mul(&lx, &rx);
    let resim = fp_mul(&ly, &ry);

    resre = fp_sub(&l0, &resre);

    (resre, resim)
}

pub fn square_and_absorb(
    fxre: &BigUint,
    fxim: &BigUint,
    ellxre: &BigUint,
    ellxim: &BigUint,
) -> (BigUint, BigUint) {
    let (resre, resim) = fp2sqr(&fxre, &fxim);
    fp2mul(&resre, &resim, &ellxre, &ellxim)
}

pub fn add_and_line(
    tx: &BigUint,
    tz: &BigUint,
    ty: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    mont_A: &BigUint,
) -> ((BigUint, BigUint, BigUint, BigUint), (BigUint, BigUint)) {
    let mut lambdax: BigUint;
    let mut lambdaz: BigUint;
    let mut ellxre: BigUint;

    lambdax = fp_mul(py, &tz);
    lambdax = fp_sub(&lambdax, &ty);

    lambdaz = fp_mul(px, &tz); //additional trick is that px can be low so implemented as adds
    lambdaz = fp_sub(&lambdaz, &tx);

    let lambdaz2 = fp_sq(&lambdaz); //tmp1 = pow(lambdaz, 2);
    let tylambdaz2 = fp_mul(&ty, &lambdaz2);

    //Get ell = line(t, p, r) = (ry - ty - lambda(rx - tx)) / (rx + tx + px - lambda*lambda + a)
    //note here that rx and px are related for specific points as sampled by us, namely rx = -px

    if *rx == fp_neg(px) {
        //something wrong here? TODO:fix
        ellxre = fp_add(&lambdaz, &tx);
        ellxre = fp_add(&ellxre, &tx);
    } else {
        ellxre = fp_mul(rx, &tz);
        ellxre = fp_sub(&tx, &ellxre);
    }

    ellxre = fp_mul(&ellxre, &lambdaz);
    ellxre = fp_mul(&ellxre, &lambdax);
    ellxre = fp_sub(&ellxre, &tylambdaz2);

    //it feels like we can improve this by using the rep (XX : XZ : ZZ : YZ)
    //right now we just use the (XZ : ZZ : YZ) part and fix it in the end with v
    //perhaps all of this can also be combined well with ellxre

    let tzlambdaz2 = fp_mul(&tz, &lambdaz2); //tmp2 = tz * tmp1;
    let lambdax2 = fp_sq(&lambdax); //tmp3 = pow(lambdax, 2);
    let txlambdaz2 = fp_mul(&tx, &lambdaz2); //tmp4 = tx * tmp1; //used in step 1.b and step 3
    let lambdaz3 = fp_mul(&lambdaz2, &lambdaz); //tmp5 = lambdaz * tmp1;
    let tzlambdax2 = fp_mul(&tz, &lambdax2); //tmp6 = tz * tmp3;

    let mut newtx = fp_add(px, mont_A);
    newtx = fp_mul(&newtx, &fp_neg(&tzlambdaz2));
    newtx = fp_sub(&newtx, &txlambdaz2);
    newtx = fp_add(&newtx, &tzlambdax2);

    let tmp = fp_mul(&lambdaz3, &ty);
    let mut newty = fp_sub(&txlambdaz2, &newtx);
    newty = fp_mul(&newty, &lambdax);
    newty = fp_sub(&newty, &tmp);

    let vx2 = fp_sq(&newtx);
    let vxz = fp_mul(&tzlambdaz2, &newtx);
    let vz2 = fp_sq(&tzlambdaz2);
    let mut vyz = fp_mul(&lambdaz, &newty);
    vyz = fp_mul(&tz, &vyz);

    ((vx2, vxz, vz2, vyz), (ellxre, tzlambdaz2))
}

pub fn give_line_bit(
    lambdax: &BigUint,
    lambday: &BigUint,
    ry: &BigUint,
) -> (BigUint, BigUint) {
    let resim = fp_mul(ry, &lambday);

    ((*lambdax).clone(), resim)
}

pub fn absorb(
    fxre: &BigUint,
    fxim: &BigUint,
    ellxre: &BigUint,
    ellxim: &BigUint,
) -> (BigUint, BigUint) {
    fp2mul(&fxre, &fxim, &ellxre, &ellxim)
}

pub fn miller(
    n: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (BigUint, BigUint) {
    //proper bigint version of the old miller4 in legacy/miller.rs
    //reworked to be highly optimized for the specific usecase

    //the point P has to be an F_p point on E and is thus given only by px and py in F_p, no F_p2 needed
    //the point R has to be an F_p point on the twist of E, and thus rx in F_p and ry is fully imaginairy, e.g. i*ry for some ry in F_p

    //returns a value in F_p^2. This value is later rased to the power (p-1) hence invariant under multiplication by F_p values
    //this also means that we can ignore the denominator when working projectively

    //this only works for the miller we are trying to compute, where p has p.x and p.y in F_p and r has p.x in F_p and p.y in i*F_p
    //in general, for CSIDH elkies, one can find a basis for E[ell_i] with such point, as both E and the twist have a point of order ell_i
    //and these generate E[ell_i]

    //this is the same as miller3, i.e. reuses computations from line for ell, is projective, and reduces to only necessary fp arith, but also uses
    //the fact that the value of ellz and fz dont matter, because they are ignored with the final raising

    //make sure we are using extraspecial form: px = -rx, ry = i*py?
    if *px == fp_neg(rx) {
        println!("Extraspecial form used");
    }
    if *py == *ry {
        println!("E0 form used");
    }

    let mut tx2: BigUint = fp_sq(px);
    let mut txtz: BigUint = px.clone();
    let mut tz2: BigUint = (*BIGONE).clone();
    let mut tytz: BigUint = py.clone();

    //our to be return value, in Fp2 but can ignore denominator, its fp invariant
    let mut fxre: BigUint = (*BIGONE).clone();
    let mut fxim: BigUint = (*BIGZERO).clone();

    //we loop per bit over n
    for i in (0..BigUint::bits(n) - 1).rev() {
        //check if t is inf or two
        if is_inf(&txtz, &tytz, &tz2) | is_two(&txtz, &tytz, &tz2) {
            //square f, set t to INF
            (fxre, fxim) = fp2sqr(&fxre, &fxim);
            set_inf(&txtz, &tytz, &tz2);
        } else {
            let lambdax: BigUint;
            let lambday: BigUint;
            let lambdaz: BigUint;

            let ellxre: BigUint;
            let ellxim: BigUint;

            //double T, get line coeffs
            ((tx2, txtz, tz2, tytz), (lambdax, lambday, lambdaz)) =
                dbl_and_line(&tx2, &txtz, &tz2, &tytz, mont_A);

            //evaluate line in rx, ry
            (ellxre, ellxim) = give_line(&lambdax, &lambday, &lambdaz, rx, ry);

            //square f and multiply by evaluation
            (fxre, fxim) = square_and_absorb(&fxre, &fxim, &ellxre, &ellxim);
        }

        //when the bit is 1
        if BigUint::bit(n, i) {
            if is_inf(&txtz, &tytz, &tz2) {
                txtz = px.clone();
                tytz = py.clone();
                tz2 = (*BIGONE).clone();
                tx2 = fp_sq(&txtz);
                continue;
            }
            let lambdax: BigUint;
            let lambday: BigUint;

            let ellxre: BigUint;
            let ellxim: BigUint;

            //if t is two, we do nothing will ell or f
            if is_two(&txtz, &tytz, &tz2) {
                //add P to T, get line coeffs, here tx2 is not actually computed
                ((tx2, txtz, tz2, tytz), (lambdax, lambday)) =
                add_and_line(&txtz, &tz2, &tytz, px, py, rx, mont_A);
                continue;
            }

            //add P to T, get line coeffs, here tx2 is not actually computed
            ((tx2, txtz, tz2, tytz), (lambdax, lambday)) =
            add_and_line(&txtz, &tz2, &tytz, px, py, rx, mont_A);

            //evaluate line in rx, ry
            (ellxre, ellxim) = give_line_bit(&lambdax, &lambday, ry);

            //just multiply f and eval
            (fxre, fxim) = absorb(&fxre, &fxim, &ellxre, &ellxim);
        }
    }

    (fxre, fxim)
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

        let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let (fxre, fxim) = miller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&fxre, &fxim, &to_red);

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

        let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        reset_big_cost();
        let (fxre, fxim) = miller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&fxre, &fxim, &to_red);
        print_big_cost();

        //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"7230A53FA464DF7559C7716823F50D114A2A96069A89F0027C01837A5BF7912642D0897BE621FADC7134A98E787C135501F462C1AEFE5DD1DDA8ED8BBE42865", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"2892CA051611D66B7A0DE4470C028535C04A57870FFA9C5FAC75644774EDE7AB2C9728E62D54CC77E33D7DC61088F4D6C873D167D0651F620C3904AA374C8C69", 16).unwrap();
        println!("{}*i + {}", zetare, zetaim);

        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }
}
