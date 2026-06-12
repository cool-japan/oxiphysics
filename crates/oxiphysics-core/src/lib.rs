// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core types, traits, and abstractions for the OxiPhysics engine.
//!
//! This crate provides fundamental building blocks used across all
//! OxiPhysics sub-crates, including vector/quaternion types, transforms,
//! bounding volumes, and core simulation traits.
//!
//! **All other OxiPhysics crates should import math types from
//! `oxiphysics_core::math` rather than depending on nalgebra directly.**
#![warn(missing_docs)]

pub mod collision;
pub mod compensated;
pub mod complex;
pub mod dual_quaternion;
pub mod exact_predicates;
pub mod extended_precision;
mod error;
pub mod graph;
pub mod interpolation;
pub mod interval;
pub mod linalg;
pub mod math;
pub mod numerical_linear_algebra;
pub mod numerics;
pub mod ode;
pub mod optimization;
pub mod parallel;
pub mod parallel_orchestrator;
pub mod random_processes;
pub mod sampling;
pub mod signal;
pub mod sparse;
pub mod spatial;
pub mod statistics;
pub mod stochastic;
pub mod tensor;
mod traits;
pub mod types;
pub mod world;

pub use dual_quaternion::*;
pub use error::*;
pub use interpolation::{
    AkimaSpline, BSplineBasis, BSplineCurve, CatmullRomSpline, Grid3Params, HermiteSpline,
    MonotoneCubicSpline, NaturalCubicSpline, NurbsCurve, RBFInterpolation, RbfKernel,
    barycentric_2d, barycentric_interp_2d, barycentric_rational, bicubic, bilinear, bilinear_grid,
    catmull_rom, dist3, hermite, hermite_deriv, hermite3, lerp, lerp3, monotone_cubic,
    natural_neighbor_interp, nurbs_evaluate_with_derivative, quat_dot, quat_nlerp, quat_norm,
    quat_normalize, quat_slerp, quat_squad, rbf_fit, rbf_interpolate, rbf_thin_plate_spline,
    rbf_tps_evaluate, rbf_tps_fit, trilinear, trilinear_grid,
};
pub use linalg::{
    characteristic_poly3, det3, frobenius_norm3, inv3, is_symmetric3, polar_decomp3, qr_decomp3,
    solve3, svd3, symmetric_eigen3, symmetric_eigenvalues3, trace3,
};
pub use math::{Mat3, Quat, Real, Vec3};
pub use numerics::*;
pub use statistics::{
    BetaDistribution, ExponentialDistribution, Histogram2D, KernelDensityEstimate,
    KernelDensityEstimate2D, MaxwellBoltzmann, NormalDistribution, PcaResult, PoissonDistribution,
    SlidingWindowStats, StatRng, UniformDistribution, WelfordOnline, acf, autocorrelation_time,
    block_average, boltzmann_factor, bootstrap_ci, bootstrap_se, chi_squared_gof,
    chi_squared_statistic, chi_squared_test, correlation, correlation_matrix, covariance,
    covariance_matrix, ecdf_at, empirical_cdf, free_energy_from_partition, histogram,
    huber_m_estimator, iqr, kruskal_wallis_h, ks_statistic, ks_test_one_sample, ks_test_two_sample,
    kurtosis, linear_regression, mad, mann_whitney_u, maxwell_boltzmann_speed, mean, median, pacf,
    partition_function, pca, pca_transform, pearson_correlation, pearson_r, percentile, quartiles,
    running_average, sample_kurtosis, sample_skewness, shapiro_wilk_w, skewness,
    spearman_correlation, std_dev, t_test_one_sample, t_test_one_sample_full, t_test_two_sample,
    tukey_biweight_estimator, variance, welch_t_test, wilcoxon_signed_rank,
};
pub use stochastic::{
    GeometricBrownianMotion, LangevinDynamics, OrnsteinUhlenbeck, RandomWalk, Rng, WienerProcess,
    autocorrelation, diffusion_coefficient, euler_maruyama, milstein,
};
pub use tensor::*;
pub use traits::*;
pub use types::{
    Aabb, BodyHandle, BoundingBox3D, Capsule, ColliderHandle, DefaultConfigs, MassProperties,
    PhysicsConfig, PhysicsConfigBuilder, Plane, Plane3D, Ray, Ray3D, RayAabbHit, Sphere, TimeStep,
    Transform, Transform3D, Triangle,
};
pub use world::{
    Body, BodyArena, BodyMode, BodyVelocity, Constraint, ConstraintKind, Island, PhysicsEvent,
    PhysicsWorld,
};
pub mod adaptive_timestepping;
pub mod autodiff;
pub mod bayesian_inference;
pub mod bayesian_opt;
pub mod cache_layout;
pub mod category_theory;
pub mod causal_inference;
pub mod chaos_theory;
pub mod complex_analysis;
pub mod compressed_sensing;
pub mod continuation;
pub mod control;
pub mod control_systems;
pub mod control_theory;
pub mod convex_analysis;
pub mod differential_geometry;
pub mod dimensionless;
pub mod dynamic_systems;
pub mod dynamical_systems;
pub mod ergodic_theory;
pub mod finite_difference;
pub mod formal_verification;
pub mod fractional_calculus;
pub mod functional_analysis;
pub mod game_theory;
pub mod geometry_primitives;
pub mod information_geometry;
pub mod information_theory;
pub mod machine_learning;
pub mod measure_theory;
pub mod metric_spaces;
pub mod monte_carlo;
pub mod multiscale_methods;
pub mod network_dynamics;
pub mod neural_ode;
pub mod numerical_integration;
pub mod numerical_methods;
pub mod numerical_methods_ext;
pub mod numerical_ode;
pub mod pde;
pub mod pde_solvers;
pub mod persistent_homology;
pub mod petri_nets;
pub mod probabilistic_models;
pub mod quadrature;
pub mod random;
pub mod renormalization_group;
pub mod signal_analysis;
pub mod signal_processing;
pub mod simd_math;
pub mod spectral_methods;
pub mod stability;
pub mod statistical_mechanics;
pub mod stochastic_processes;
pub mod symbolic_algebra;
pub mod symbolic_math;
pub mod time_series_analysis;
pub mod topological_data;
pub mod topology;
pub mod topology_optimization;
pub mod uncertainty_quantification;
pub mod unit_conversion;
pub mod wavelet_transform;
pub mod wavelets;

#[cfg(feature = "scirs2")]
pub mod scirs2_integrator;
