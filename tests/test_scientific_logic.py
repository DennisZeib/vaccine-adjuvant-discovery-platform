import numpy as np
import pandas as pd

from app import (
    compute_regression_metrics,
    fit_calibration_model,
    pareto_front,
    run_dfba_sanity_tests,
)


def test_compute_regression_metrics_returns_expected_keys():
    metrics = compute_regression_metrics(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.1]))
    assert set(metrics.keys()) >= {"rmse", "r2"}
    assert metrics["rmse"] < 0.2
    assert metrics["r2"] > 0.95


def test_fit_calibration_model_returns_train_and_test_metrics():
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    y = np.array([0.2, 0.35, 0.55, 0.85, 1.3, 2.0, 3.0, 4.6])

    result = fit_calibration_model(t, y, "Exponential", test_size=0.25, random_state=7)

    assert result["params"]["model"] == "exponential"
    assert result["params"]["mu"] > 0
    assert len(result["train_idx"]) > 0
    assert len(result["test_idx"]) > 0
    assert set(result["train_metrics"].keys()) >= {"rmse", "r2"}
    assert set(result["test_metrics"].keys()) >= {"rmse", "r2"}


def test_pareto_front_keeps_only_non_dominated_points():
    df = pd.DataFrame({
        "avg_score": [90.0, 80.0, 70.0],
        "toxicity": [1.0, 2.0, 3.0],
    })

    mask = pareto_front(df, "avg_score", "toxicity")

    assert mask.tolist() == [True, False, False]


def test_run_dfba_sanity_tests_returns_expected_structure():
    result = run_dfba_sanity_tests(
        S0_mM=1.0,
        duration_h=2.0,
        dt_h=0.5,
        Vmax_a=10.0,
        Km_a=0.5,
        Vmax_b=10.0,
        Km_b=0.5,
        biomass_conv=0.001,
        volume_l=1.0,
    )

    assert set(result.keys()) >= {"unlimited", "limited"}
    assert "ok" in result["unlimited"]
    assert "ok" in result["limited"]
