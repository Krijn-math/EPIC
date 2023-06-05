//using ideas from 2019/077 we separate the calculation of the line function and the evaluation of the pairing in a specific point
//this requires a bit of memory, namely log(p)^2 but is worth it in the search of full torsion points
//in honor of Michael Scott, I've named it the Scott-separated Miller algorithm

//I did so by declaring a line structure, and then an array of such lines
//this array of lines is computed a la Algorithm 2 of 2019/077
//and can be evaluated a la Algorithm 3 of 2019/077
//mathematically, the array of lines expresses f_{rP}(-) = e(P, -) as the miller function E^t --> mu_p+1
//and we simply evaluate these in e(P, Q) in algorithm 3

use crate::bigconstants::*;
use crate::bigint::*;
use num_bigint::BigUint;
use super::addchain::ADDCHAIN;
use super::evaluationhelperfunctions::*;
use super::scott_window::*;
use super::millerfunctions::*;
use super::constructionhelperfunctions::*;

pub fn scottmiller_construction(
    addchain: &ADDCHAIN,
    px: &BigUint,
    py: &BigUint,
    ppre: Vec<AffPoint>,
    mont_A: &BigUint,
) -> Vec<Line> {
    let mut t: Point = point_cast(&px, &py);
    let mut f: Vec<Line> = Default::default();
    let mut step: Line;

    for i in 0..addchain.len() {
        let j = addchain[i];
        if j == 0 {
            (t, step) = dbl_and_line(t, mont_A);
            f.push(step);
            break;
        } 
        
        (t, step) = add_and_line(t, &ppre[j], mont_A);              //can be substraction, depends on j
        f.push(step);
    }

    f
}

pub fn scottmiller_evaluation(
    addchain: &ADDCHAIN,
    f: Vec<Line>,
    respre: Vec<Result>,
    rx: &BigUint,
    ry: &BigUint,
)   -> Result {
    let mut re: BigUint = (*BIGONE).clone();
    let mut im: BigUint = (*BIGZERO).clone();
    let mut res = Result { re, im };

    //here should be the evaluation loop 
    for i in 0..addchain.len() {
        let j = addchain[i];
        if j == 0 {
            res = square(res);
            let eval = apply(f[i].clone(), &rx, &ry);            //step[i] = ell_{T, T}(-) so evaluate in R
            res = absorb(res, eval);
        } else {                                         //positive means add, negative means sub
            let eval = apply(f[i].clone(), &rx, &ry);             //step[i] = ell_{T, jP}(-) so evaluate in R
            res = absorb(res, eval);                                 //res <-- res*g[i](Q)
            res = absorb(res, respre[j].clone());                       //multiply by f_{jP}(Q) which is precomputed
        }
    }
    res
}

pub fn scottmiller(
    addchain: &ADDCHAIN,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
)   -> (BigUint, BigUint)
 {
    //example of how a single pairing would work
    let (ppre, fpre) = precompute_const(px, py, mont_A);
    let f = scottmiller_construction(&addchain, px, py, ppre, mont_A);

    let respre = precompute_eval(fpre, rx, ry);
    let res = scottmiller_evaluation(&addchain, f, respre, rx, ry);

    (res.re, res.im)
}


#[cfg(test)]
mod tests {
    use crate::scottmiller::addchain::ADDCHAIN;

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

        let (resre, resim) = scottmiller(&ADDCHAIN, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&resre, &resim, &to_red);

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

        let (resre, resim) = scottmiller(&ADDCHAIN, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (zetare, zetaim) = fp2pow(&resre, &resim, &to_red);

        //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"7230A53FA464DF7559C7716823F50D114A2A96069A89F0027C01837A5BF7912642D0897BE621FADC7134A98E787C135501F462C1AEFE5DD1DDA8ED8BBE42865", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"2892CA051611D66B7A0DE4470C028535C04A57870FFA9C5FAC75644774EDE7AB2C9728E62D54CC77E33D7DC61088F4D6C873D167D0651F620C3904AA374C8C69", 16).unwrap();

        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }
}
