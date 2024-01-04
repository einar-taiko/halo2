use std::fmt::Debug;

use super::commitment::{KZGCommitmentScheme, ParamsKZG};
use crate::{
    arithmetic::{parallelize, CurveAffine},
    poly::commitment::MSM,
    ZalRef,
};
use group::{Curve, Group};
use halo2curves::{
    pairing::{Engine, MillerLoopResult, MultiMillerLoop},
    zal::{self, H2cEngine, MsmAccel},
};

/// A multiscalar multiplication in the polynomial commitment scheme
#[derive(Clone, Debug)]
pub struct MSMKZG<'zal, E: Engine> {
    pub(crate) zal: ZalRef<'zal>,
    pub(crate) scalars: Vec<E::Scalar>,
    pub(crate) bases: Vec<E::G1>,
}

impl<'zal, E: Engine> MSMKZG<'zal, E> {
    /// Create an empty MSM instance
    pub fn new(zal: ZalRef<'zal>) -> Self {
        MSMKZG {
            zal,
            scalars: vec![],
            bases: vec![],
        }
    }

    /// Prepares all scalars in the MSM to linear combination
    pub fn combine_with_base(&mut self, base: E::Scalar) {
        use ff::Field;
        let mut acc = E::Scalar::ONE;
        if !self.scalars.is_empty() {
            for scalar in self.scalars.iter_mut().rev() {
                *scalar *= &acc;
                acc *= base;
            }
        }
    }
}

impl<'zal, E: Engine + Debug> MSM<E::G1Affine> for MSMKZG<'zal,E> {
    fn append_term(&mut self, scalar: E::Scalar, point: E::G1) {
        self.scalars.push(scalar);
        self.bases.push(point);
    }

    fn add_msm(&mut self, other: &Self) {
        self.scalars.extend(other.scalars().iter());
        self.bases.extend(other.bases().iter());
    }

    fn scale(&mut self, factor: E::Scalar) {
        if !self.scalars.is_empty() {
            parallelize(&mut self.scalars, |scalars, _| {
                for other_scalar in scalars {
                    *other_scalar *= &factor;
                }
            })
        }
    }

    fn check(&self) -> bool {
        bool::from(self.eval().is_identity())
    }

    fn eval(&self) -> E::G1 {
        use group::prime::PrimeCurveAffine;
        let mut bases = vec![E::G1Affine::identity(); self.scalars.len()];
        E::G1::batch_normalize(&self.bases, &mut bases);
        let engine = H2cEngine::new();
        engine.msm(&self.scalars, &bases)
    }

    fn bases(&self) -> Vec<E::G1> {
        self.bases.clone()
    }

    fn scalars(&self) -> Vec<E::Scalar> {
        self.scalars.clone()
    }
}

/// A projective point collector
#[derive(Debug, Clone)]
pub(crate) struct PreMSM<'zal, E: Engine> {
    projectives_msms: Vec<MSMKZG<'zal, E>>,
    zal: ZalRef<'zal>,
}

impl<'zal, E: Engine + Debug> PreMSM<'zal, E> {
    pub(crate) fn new(zal: ZalRef<'zal>) -> Self {
        PreMSM {
            projectives_msms: vec![],
            zal,
        }
    }

    pub(crate) fn normalize(self) -> MSMKZG<'zal, E> {
        use group::prime::PrimeCurveAffine;

        let (scalars, bases) = self
            .projectives_msms
            .into_iter()
            .map(|msm| (msm.scalars, msm.bases))
            .unzip::<_, _, Vec<_>, Vec<_>>();

        MSMKZG {
            zal: self.zal,
            scalars: scalars.into_iter().flatten().collect(),
            bases: bases.into_iter().flatten().collect(),
        }
    }

    pub(crate) fn add_msm(&mut self, other: MSMKZG<E>) {
        self.projectives_msms.push(other);
    }
}

/*
impl<'params, 'zal, E: MultiMillerLoop + Debug> From<&'params ParamsKZG<E>>
    for DualMSM<'params, 'zal, E>
{
    fn from(params: &'params ParamsKZG<E>, zal: ZalRef<'zal>) -> Self {
        DualMSM::new(params, zal)
    }
}
*/

/// Two channel MSM accumulator
#[derive(Debug, Clone)]
pub struct DualMSM<'a, 'zal, E: Engine> {
    pub(crate) params: &'a ParamsKZG<E>,
    pub(crate) left: MSMKZG<'zal, E>,
    pub(crate) right: MSMKZG<'zal, E>,
}

impl<'a, 'zal, E: MultiMillerLoop + Debug> DualMSM<'a, 'zal, E> {
    /// Create a new two channel MSM accumulator instance
    pub fn new(params: &'a ParamsKZG<E>, zal: ZalRef<'zal>) -> Self {
        Self {
            params,
            left: MSMKZG::new(zal),
            right: MSMKZG::new(zal),
        }
    }

    /// Scale all scalars in the MSM by some scaling factor
    pub fn scale(&mut self, e: E::Scalar) {
        self.left.scale(e);
        self.right.scale(e);
    }

    /// Add another multiexp into this one
    pub fn add_msm(&mut self, other: Self) {
        self.left.add_msm(&other.left);
        self.right.add_msm(&other.right);
    }

    /// Performs final pairing check with given verifier params and two channel linear combination
    pub fn check(self) -> bool {
        let s_g2_prepared = E::G2Prepared::from(self.params.s_g2);
        let n_g2_prepared = E::G2Prepared::from(-self.params.g2);

        let left = self.left.eval();
        let right = self.right.eval();

        let (term_1, term_2) = (
            (&left.into(), &s_g2_prepared),
            (&right.into(), &n_g2_prepared),
        );
        let terms = &[term_1, term_2];

        bool::from(
            E::multi_miller_loop(&terms[..])
                .final_exponentiation()
                .is_identity(),
        )
    }
}
