# Global Alpha Model for Interest Rate Spread Options

## 1. Objective

In the interest rate (IR) spread options market, multiple underlying rates are traded, each with their own volatility term structures and smiles.  
The goal is to jointly model these underlyings in a **consistent multi-factor stochastic volatility framework** that:
- Captures realistic smile shapes across all markets.
- Ensures **no-arbitrage** in both strike and maturity.
- Embeds cross-asset correlation naturally into the dynamics.

---

## 2. Multi-Asset Core Structure

We extend a **Heston-type stochastic volatility process** to multiple underlyings:

\[
dS_i = \sum_{j=1}^m \alpha_{ij} v_j S_i \, dW_j + v_{m+i} S_i \, dW_{m+i}
\]
\[
dv_j = \lambda_j (\theta_j - v_j) dt + \nu_j v_j \, dW_{noise(j)}
\]

Where:
- **\( S_i \)** = underlying rate (e.g., swap rate for a given tenor)
- **\( v_j \)** = volatility factors (some shared across markets, others idiosyncratic)
- **\( \alpha_{ij} \)** = factor loadings linking vol factors to underlyings
- **\( \rho \)** = instantaneous correlations between returns and volatility shocks

This shared factor structure ensures coherent dynamics across all underlyings.

---

## 3. Arbitrage-Free Smile Parameterization

For each underlying, the model produces three key **observables**:
1. **V0** — Instantaneous variance.
2. **VarVar** — Variance of volatility.
3. **CoVar** — Covariance between underlyings' volatilities.

The smile is parameterized directly from these observables to guarantee:
- **No calendar arbitrage** (consistency across maturities).
- **No butterfly arbitrage** (convexity in strike space).
- Smooth, stable implied volatility surfaces.

Because these parameters come from the dynamics, the smile is **not ad hoc** but economically grounded.

---

## 4. Spread Option Pricing

Pricing steps:
1. **Compute observables** for each underlying.
2. **Correct for \( \nu \to 0 \)** limit to ensure proper small-vol-of-vol behavior.
3. Apply **Normal Heston** framework in multi-asset form to obtain spread option prices.

Correlations between underlyings are given by the factor structure, ensuring **joint arbitrage-free pricing**.

---

## 5. Benefits in IR Spread Markets

- **Unified calibration**: All related IR markets calibrated together.
- **Arbitrage-free** smiles across strikes and maturities.
- **Robust** from short to long maturities.
- **Transparent**: Parameters have direct stochastic interpretations.

---

## 6. Conclusion

The Global Alpha Model provides a **transparent, modular, and consistent** framework for modeling and pricing IR spread options with realistic and arbitrage-free volatility smiles.
