/* ============================================================================
 * truth-recovery/dgp-ipd.mjs
 * STANDALONE seeded known-truth DGP for binary IPD meta-analysis.
 *
 * Truth model (random-effects logistic, the canonical IPD-MA generative model):
 *   For study i = 1..k:
 *     alpha_i  = alpha0 + u_i,   u_i ~ N(0, sigmaAlpha^2)   (baseline-risk variation)
 *     theta_i  = theta  + v_i,   v_i ~ N(0, tau^2)          (treatment-effect heterogeneity)
 *   For patient j in study i, arm t in {0,1}:
 *     logit P(Y=1) = alpha_i + theta_i * t
 *
 * KNOWN TRUTH = overall treatment log-OR `theta` and between-study SD `tau`.
 * Individual binary outcomes are produced. Output rows: {study, treatment, y}
 * -- exactly the shape OneStageGLMM.fitBinary and a 2x2-per-study two-stage
 * pooling both consume.
 * ==========================================================================*/

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function gauss(rng) {
  let u1 = rng(); while (u1 <= 0) u1 = rng();
  const u2 = rng();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}
function sigmoid(x) { return 1 / (1 + Math.exp(-x)); }

/**
 * Generate one known-truth IPD dataset.
 * @returns {rows:[{study,treatment,y}], truth:{theta,tau,k,nPerArm}}
 */
function generateIPD({ seed, k, nPerArm, alpha0, theta, sigmaAlpha, tau }) {
  const rng = mulberry32(seed);
  const rows = [];
  for (let i = 0; i < k; i++) {
    const alpha_i = alpha0 + sigmaAlpha * gauss(rng);
    const theta_i = theta + tau * gauss(rng);
    const sid = 'S' + (i + 1);
    for (let t = 0; t <= 1; t++) {
      const p = sigmoid(alpha_i + theta_i * t);
      for (let j = 0; j < nPerArm; j++) {
        const y = rng() < p ? 1 : 0;
        rows.push({ study: sid, treatment: t, y: y });
      }
    }
  }
  return { rows, truth: { theta, tau, k, nPerArm, alpha0, sigmaAlpha } };
}

/**
 * Per-study 2x2 -> log-OR + variance (Woolf, with 0.5 continuity correction
 * only when a zero cell is present -- per advanced-stats rule). For two-stage.
 * @returns {effects:[logOR_i], variances:[v_i], tables:[{a,b,c,d}]}
 */
function perStudyLogOR(rows) {
  const byStudy = new Map();
  for (const r of rows) {
    if (!byStudy.has(r.study)) byStudy.set(r.study, { a: 0, b: 0, c: 0, d: 0 });
    const t = byStudy.get(r.study);
    // a=treated event, b=treated no-event, c=control event, d=control no-event
    if (r.treatment === 1) { if (r.y === 1) t.a++; else t.b++; }
    else { if (r.y === 1) t.c++; else t.d++; }
  }
  const effects = [], variances = [], tables = [];
  for (const t of byStudy.values()) {
    let { a, b, c, d } = t;
    if (a === 0 || b === 0 || c === 0 || d === 0) { a += 0.5; b += 0.5; c += 0.5; d += 0.5; }
    const logOR = Math.log((a * d) / (b * c));
    const v = 1 / a + 1 / b + 1 / c + 1 / d;
    effects.push(logOR); variances.push(v); tables.push({ a, b, c, d });
  }
  return { effects, variances, tables };
}

export { mulberry32, gauss, sigmoid, generateIPD, perStudyLogOR };
