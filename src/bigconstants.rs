use lazy_static::lazy_static;
use num_bigint::{BigUint, ToBigUint, BigInt, ToBigInt};
use num_traits::{One, Zero}; // 1.4.0

lazy_static! {
     //CSIDH-44
    // pub static ref CSIDH512: BigUint = BigUint::parse_bytes(b"14841476269619", 10).unwrap();  //4*prod of [ 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 ]
    // pub static ref ELLS: Vec<i128> = [ 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 ].to_vec();

    //CSIDH512
    pub static ref CSIDH512: BigUint = BigUint::parse_bytes(b"65B48E8F740F89BFFC8AB0D15E3E4C4AB42D083AEDC88C425AFBFCC69322C9CDA7AAC6C567F35507516730CC1F0B4F25C2721BF457ACA8351B81B90533C6C87B", 16).unwrap();
    pub static ref ELLS: Vec<i128> = [ 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                                            73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157,
                                            163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241,
                                            251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
                                            349, 353, 359, 367, 373, 587 ].to_vec();

    pub static ref NELLS: usize = ELLS.len();

    pub static ref BIGZERO: BigUint = BigUint::zero();
    pub static ref BIGONE: BigUint =  BigUint::one();
    pub static ref BIGTWO: BigUint = 2.to_biguint().unwrap();

    pub static ref INV: BigUint = &*CSIDH512 - &(*BIGTWO);
    pub static ref PPLUSONE: BigUint = (&*CSIDH512) + &(*BIGONE);
    pub static ref LEGENDRE: BigUint = (&*CSIDH512 - &(*BIGONE)) >> 1;
    pub static ref SQRT: BigUint = (&*CSIDH512 + &(*BIGONE)) >> 2;
    pub static ref TRICK: BigUint = &*CSIDH512 + &*CSIDH512 + &*BIGTWO;

    pub static ref BOUND: BigUint = BigUint::parse_bytes(b"2856F1399D91D6592142B95420000000000000000000000000000000000000000", 16).unwrap();

    pub static ref COST_DBL_MUL: usize = 3;
    pub static ref COST_DBL_SQR: usize = 2;
    pub static ref COST_ADD_MUL: usize = 4;
    pub static ref COST_ADD_SQR: usize = 2;
    pub static ref COST_LDD_MUL: usize = *COST_ADD_MUL + *COST_DBL_MUL;
    pub static ref COST_LDD_SQR: usize = *COST_ADD_SQR + *COST_DBL_SQR;
    pub static ref COST_SPE_MUL: usize = *COST_ADD_MUL + *COST_DBL_MUL - 1;
    pub static ref COST_SPE_SQR: usize = *COST_ADD_SQR + *COST_DBL_SQR;

    pub static ref NAF: Vec<i8> = naf(&(*PPLUSONE).clone().to_bigint().unwrap());
}

pub fn log2(n: usize) -> usize {
    (n.count_ones() + n.count_zeros() - n.leading_zeros()) as usize
}

pub fn cofac() -> usize {
    //assumes cofac is a power of two
    let mut z: BigUint = &(*CSIDH512) + &(*BIGONE);
    let mut res: usize = 0;
    while &z % &*BIGTWO == *BIGZERO {
        res += 1;
        z /= &*BIGTWO;
    }
    res
}

pub fn ellbits() -> usize {
    let z: BigUint = (&(*CSIDH512) + &(*BIGONE)) >> cofac();
    BigUint::bits(&z) as usize
}

pub fn is_zero(n: &BigUint) -> bool {
    BigUint::is_zero(n)
}

pub fn is_one(n: &BigUint) -> bool {
    BigUint::is_one(n)
}

pub static mut BIGCOST: [u128; 3] = [0, 0, 0];

pub fn reset_big_cost() {
    unsafe {
        BIGCOST = [0, 0, 0];
    }
}

pub fn print_big_cost() {
    unsafe {
        println!(
            "MULLS {}, SQUARES {}, ADDS {}",
            BIGCOST[0], BIGCOST[1], BIGCOST[2]
        );
    }
}

pub fn mid(upp: usize, low: usize) -> usize {
    if (upp + low) % 2 == 0 {
        return (upp + low) / 2;
    }

    (upp + low - 1) / 2
}

pub fn ProdVec(v: &[i128], bot: usize, top: usize) -> BigUint {
    let mut res: BigUint = BigUint::one();
    let mut i = bot;
    while i <= top {
        res *= v[i].to_biguint().unwrap();
        i += 1;
    }

    res
}

pub fn naf(y: &BigInt) -> Vec<i8> {
    //computes naf, given by  Prodinger, Helmut. "On Binary Representations of Integers with Digits -1, 0, 1"

    let x = (*y).clone();
    let xh = &x >> 1;
    let x3 = &x + &xh;
    let c = &xh ^ &x3;
    let np = &x3 & &c;
    let nm = &xh & &c;

    let zeroes = (BigInt::bits(&np) - BigInt::bits(&nm)) as usize;
    let mut vec = vec![0; zeroes];

    let vecnp = np.to_radix_be(2).1;
    let mut vecnm = nm.to_radix_be(2).1;
    vec.append(&mut vecnm);

    let ln = BigInt::bits(&np) as usize;

    let mut vec2: Vec<i8> = vec![0; ln];

    for i in 0..ln {
        vec2[i] = (vecnp[i] as i8) - (vec[i] as i8);
        //print!("{}", (vecnp[i] as i8) + 2*(vec[i] as i8));
    }

    vec2
}
