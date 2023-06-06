//here, I've tried to implement the approach from 2016/963 to optimize the arith for T <- T + T and (fxre, fxim)

use crate::bigconstants::*;
use crate::bigint::*;
use crate::tradmiller::bigmiller::*;
use num_bigint::BigUint;
use num_bigint::ToBigInt;

pub fn double(
    tx2: &BigUint,
    txtz: &BigUint,
    tz2: &BigUint,
    tytz: &BigUint,
    fxre: &BigUint,
    fxim: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> ((BigUint, BigUint, BigUint, BigUint), (BigUint, BigUint)) {
    //check if t is inf or two
    if is_inf(txtz, tytz, tz2) | is_two(txtz, tytz, tz2) {
        //square f, set t to INF
        let (resfxre, resfxim) = fp2sqr(fxre, fxim);
        let restx2 = (*tx2).clone();
        let (restxtz, restytz, restz2) =
            ((*BIGZERO).clone(), (*BIGONE).clone(), (*BIGZERO).clone());
        ((restx2, restxtz, restz2, restytz), (resfxre, resfxim))
    } else {
        //double T, get line coeffs
        let ((restx2, restxtz, restz2, restytz), (lambdax, lambday, lambdaz)) =
            dbl_and_line(tx2, txtz, tz2, tytz, mont_A);

        //evaluate line in rx, ry
        let (ellxre, ellxim) = give_line(&lambdax, &lambday, &lambdaz, rx, ry);

        //square f and multiply by evaluation
        let (resfxre, resfxim) = square_and_absorb(fxre, fxim, &ellxre, &ellxim);
        ((restx2, restxtz, restz2, restytz), (resfxre, resfxim))
    }
}

pub fn add(
    _tx2: &BigUint,
    txtz: &BigUint,
    tz2: &BigUint,
    tytz: &BigUint,
    fxre: &BigUint,
    fxim: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> ((BigUint, BigUint, BigUint, BigUint), (BigUint, BigUint)) {
    let restx2: BigUint;
    let restxtz: BigUint;
    let restz2: BigUint;
    let restytz: BigUint;

    if is_inf(txtz, tytz, tz2) {
        restxtz = px.clone();
        restytz = py.clone();
        restz2 = (*BIGONE).clone();
        restx2 = fp_sq(txtz);
        let resfxre = (*fxre).clone();
        let resfxim = (*fxim).clone();
        return ((restx2, restxtz, restz2, restytz), (resfxre, resfxim));
    }

    //if t is two, we do nothing will ell or f
    if is_two(txtz, tytz, tz2) {
        //add P to T, get line coeffs, here tx2 is not actually computed
        let ((restx2, restxtz, restz2, restytz), (_lambdax, _lambday)) =
            add_and_line(txtz, tz2, tytz, px, py, rx, mont_A);
        let resfxre = (*fxre).clone();
        let resfxim = (*fxim).clone();
        return ((restx2, restxtz, restz2, restytz), (resfxre, resfxim));
    }

    //add P to T, get line coeffs, here tx2 is not actually computed
    let ((restx2, restxtz, restz2, restytz), (lambdax, lambday)) =
        add_and_line(txtz, tz2, tytz, px, py, rx, mont_A);
    let (ellxre, ellxim) = give_line_bit(&lambdax, &lambday, ry);
    let (resfxre, resfxim) = absorb(fxre, fxim, &ellxre, &ellxim);

    ((restx2, restxtz, restz2, restytz), (resfxre, resfxim))
}

pub fn subtract(
    tx2: &BigUint,
    txtz: &BigUint,
    tz2: &BigUint,
    tytz: &BigUint,
    fxre: &BigUint,
    fxim: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> ((BigUint, BigUint, BigUint, BigUint), (BigUint, BigUint)) {
    add(
        tx2, txtz, tz2, tytz, fxre, fxim, px, &fp_neg(py), rx, ry, mont_A,
    )
}


pub fn cleanmiller(
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

    let mut tx2: BigUint = fp_sq(px);
    let mut txtz: BigUint = px.clone();
    let mut tz2: BigUint = (*BIGONE).clone();
    let mut tytz: BigUint = py.clone();

    //our to be return value, in Fp2 but can ignore denominator, its fp invariant
    let mut fxre: BigUint = (*BIGONE).clone();
    let mut fxim: BigUint = (*BIGZERO).clone();

    //we loop per bit over n
    for i in (0..BigUint::bits(n) - 1).rev() {
        ((tx2, txtz, tz2, tytz), (fxre, fxim)) =
            double(&tx2, &txtz, &tz2, &tytz, &fxre, &fxim, rx, ry, mont_A);

        //when the bit is 1
        if BigUint::bit(n, i) {
            ((tx2, txtz, tz2, tytz), (fxre, fxim)) = add(
                &tx2, &txtz, &tz2, &tytz, &fxre, &fxim, px, py, rx, ry, mont_A,
            );
        }
    }



    (fxre, fxim)
}

pub fn nafmiller(
    n: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (BigUint, BigUint) {
    //cleanmiller but implemented using a NAF, using the fact that we can easily add or substract

    let mut tx2: BigUint = fp_sq(px);
    let mut txtz: BigUint = px.clone();
    let mut tz2: BigUint = (*BIGONE).clone();
    let mut tytz: BigUint = py.clone();

    let nafn = naf(&(*n).to_bigint().unwrap());

    //our to be return value, in Fp2 but can ignore denominator, its fp invariant
    let mut fxre: BigUint = (*BIGONE).clone();
    let mut fxim: BigUint = (*BIGZERO).clone();

    //we loop per bit over n
    for _i in 1..nafn.len() {
        ((tx2, txtz, tz2, tytz), (fxre, fxim)) =
            double(&tx2, &txtz, &tz2, &tytz, &fxre, &fxim, rx, ry, mont_A);

        //add or subtract depending on the bit
        if nafn[_i] == 1 {
            ((tx2, txtz, tz2, tytz), (fxre, fxim)) = add(
                &tx2, &txtz, &tz2, &tytz, &fxre, &fxim, px, py, rx, ry, mont_A,
            );
        } else if nafn[_i] == -1 {
            ((tx2, txtz, tz2, tytz), (fxre, fxim)) = subtract(
                &tx2, &txtz, &tz2, &tytz, &fxre, &fxim, px, py, rx, ry, mont_A,
            );
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

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let (fxre, fxim) = nafmiller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
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

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let (fxre, fxim) = nafmiller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&fxre, &fxim, &to_red);

        //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"7230A53FA464DF7559C7716823F50D114A2A96069A89F0027C01837A5BF7912642D0897BE621FADC7134A98E787C135501F462C1AEFE5DD1DDA8ED8BBE42865", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"2892CA051611D66B7A0DE4470C028535C04A57870FFA9C5FAC75644774EDE7AB2C9728E62D54CC77E33D7DC61088F4D6C873D167D0651F620C3904AA374C8C69", 16).unwrap();

        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }
}
