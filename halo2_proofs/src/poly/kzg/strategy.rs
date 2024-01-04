use std::{fmt::Debug, marker::PhantomData};

use super::{
    commitment::{KZGCommitmentScheme, ParamsKZG},
    msm::{DualMSM, MSMKZG},
    multiopen::VerifierGWC,
};
use crate::{
    helpers::SerdeCurveAffine,
    plonk::Error,
    poly::{
        commitment::{Verifier, MSM},
        ipa::msm::MSMIPA,
        strategy::{Guard, VerificationStrategy},
    },
    transcript::{EncodedChallenge, TranscriptRead},
    ZalRef,
};
use ff::{Field, PrimeField};
use group::Group;
use halo2curves::{
    pairing::{Engine, MillerLoopResult, MultiMillerLoop},
    CurveAffine,
};
use rand_core::OsRng;

/// Wrapper for linear verification accumulator
#[derive(Debug, Clone)]
pub struct GuardKZG<'params, 'zal, E: MultiMillerLoop + Debug> {
    pub(crate) msm_accumulator: DualMSM<'params, 'zal, E>,
}

/// Define accumulator type as `DualMSM`
impl<'params, 'zal, E> Guard<KZGCommitmentScheme<E>> for GuardKZG<'params, 'zal, E>
where
    E::Scalar: PrimeField,
    E: MultiMillerLoop + Debug,
    E::G1Affine: SerdeCurveAffine,
    E::G2Affine: SerdeCurveAffine,
{
    type MSMAccumulator = DualMSM<'params, 'zal, E>;
}

/// KZG specific operations
impl<'params, 'zal, E: MultiMillerLoop + Debug> GuardKZG<'params, 'zal, E> {
    pub(crate) fn new(msm_accumulator: DualMSM<'params, 'zal, E>) -> Self {
        Self { msm_accumulator }
    }
}

/// A verifier that checks multiple proofs in a batch
#[derive(Clone, Debug)]
pub struct AccumulatorStrategy<'params, 'zal, E: Engine> {
    pub(crate) msm_accumulator: DualMSM<'params, 'zal, E>,
}

impl<'params, 'zal, E: MultiMillerLoop + Debug> AccumulatorStrategy<'params, 'zal, E> {
    /// Constructs an empty batch verifier
    pub fn new(params: &'params ParamsKZG<E>, zal: ZalRef) -> Self {
        AccumulatorStrategy {
            msm_accumulator: DualMSM::new(params, zal),
        }
    }

    /// Constructs and initialized new batch verifier
    pub fn with(msm_accumulator: DualMSM<'params, 'zal, E>) -> Self {
        AccumulatorStrategy { msm_accumulator }
    }
}

/// A verifier that checks a single proof
#[derive(Clone, Debug)]
pub struct SingleStrategy<'params, 'zal, E: Engine> {
    pub(crate) msm: DualMSM<'params, 'zal, E>,
}

impl<'params, 'zal, E: MultiMillerLoop + Debug> SingleStrategy<'params, 'zal, E> {
    /// Constructs an empty batch verifier
    pub fn new(params: &'params ParamsKZG<E>, zal: ZalRef) -> Self {
        SingleStrategy {
            msm: DualMSM::new(params, zal),
        }
    }
}

impl<
        'params,
        'zal,
        E: MultiMillerLoop + Debug,
        V: Verifier<
            'params,
            'zal,
            KZGCommitmentScheme<E>,
            MSMAccumulator = DualMSM<'params, 'zal, E>,
            Guard = GuardKZG<'params, 'zal, E>,
        >,
    > VerificationStrategy<'params, 'zal, KZGCommitmentScheme<E>, V>
    for AccumulatorStrategy<'params, 'zal, E>
where
    E::Scalar: PrimeField,
    E::G1Affine: SerdeCurveAffine,
    E::G2Affine: SerdeCurveAffine,
{
    type Output = Self;

    fn new(params: &'params ParamsKZG<E>, zal: ZalRef) -> Self {
        AccumulatorStrategy::new(params, zal)
    }

    fn process(
        mut self,
        f: impl FnOnce(V::MSMAccumulator) -> Result<V::Guard, Error>,
    ) -> Result<Self::Output, Error> {
        self.msm_accumulator.scale(E::Scalar::random(OsRng));

        // Guard is updated with new msm contributions
        let guard = f(self.msm_accumulator)?;
        Ok(Self {
            msm_accumulator: guard.msm_accumulator,
        })
    }

    fn finalize(self) -> bool {
        self.msm_accumulator.check()
    }
}

impl<
        'params,
        'zal,
        E: MultiMillerLoop + Debug,
        V: Verifier<
            'params,
            'zal,
            KZGCommitmentScheme<E>,
            MSMAccumulator = DualMSM<'params, 'zal, E>,
            Guard = GuardKZG<'params, 'zal, E>,
        >,
    > VerificationStrategy<'params, 'zal, KZGCommitmentScheme<E>, V>
    for SingleStrategy<'params, 'zal, E>
where
    E::Scalar: PrimeField,
    E::G1Affine: SerdeCurveAffine,
    E::G2Affine: SerdeCurveAffine,
{
    type Output = ();

    fn new(params: &'params ParamsKZG<E>, zal: ZalRef) -> Self {
        Self::new(params, zal)
    }

    fn process(
        self,
        f: impl FnOnce(V::MSMAccumulator) -> Result<V::Guard, Error>,
    ) -> Result<Self::Output, Error> {
        // Guard is updated with new msm contributions
        let guard = f(self.msm)?;
        let msm = guard.msm_accumulator;
        if msm.check() {
            Ok(())
        } else {
            Err(Error::ConstraintSystemFailure)
        }
    }

    fn finalize(self) -> bool {
        unreachable!();
    }
}
