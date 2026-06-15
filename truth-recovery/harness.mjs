/* ============================================================================
 * truth-recovery/harness.mjs
 * Wires the standalone known-truth binary IPD DGP into the REPO'S OWN
 * one-stage GLMM (OneStageGLMM.fitBinary, PQL with random study intercepts)
 * AND two-stage pooling (MetaAnalysis.randomEffectsREML + applyHKSJ on
 * per-study log-ORs), then measures bias and coverage of the known true log-OR.
 *
 * TRUTH-FIRST: every number is measured from the repo's own code on
 * known-truth inputs. Nothing hand-tuned.
 * ==========================================================================*/

import { MetaAnalysis, OneStageGLMM } from './engine.mjs';
import { generateIPD, perStudyLogOR } from './dgp-ipd.mjs';

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }

/** One replicate: returns {oneStage:{est,se,lo,hi}, twoStageZ:{...}, twoStageHKSJ:{...}} or null on failure. */
function fitOnce(rows) {
  const out = {};

  // ---- ONE-STAGE: repo's OWN GLMM (random study intercept via PQL) ----
  const os = OneStageGLMM.fitBinary(rows, 'y', 'treatment', 'study', [], {});
  if (os && !os.error && isFinite(os.pooled_effect) && isFinite(os.SE) && os.SE > 0) {
    out.oneStage = {
      est: os.pooled_effect, se: os.SE,
      lo: os.pooled_effect - 1.96 * os.SE, hi: os.pooled_effect + 1.96 * os.SE,
    };
  }

  // ---- TWO-STAGE: per-study log-OR -> repo's OWN REML RE pooling ----
  const { effects, variances } = perStudyLogOR(rows);
  const re = MetaAnalysis.randomEffectsREML(effects, variances);
  if (re && isFinite(re.pooled) && isFinite(re.se) && re.se > 0) {
    // naive z CI
    out.twoStageZ = {
      est: re.pooled, se: re.se,
      lo: re.pooled - 1.96 * re.se, hi: re.pooled + 1.96 * re.se,
    };
    // HKSJ-adjusted CI (repo's OWN applyHKSJ -> uses t_{k-1})
    const hk = MetaAnalysis.applyHKSJ(re, effects, variances);
    const lo = isFinite(hk.lowerHKSJ) ? hk.lowerHKSJ : (hk.ciHKSJ && hk.ciHKSJ[0]);
    const hi = isFinite(hk.upperHKSJ) ? hk.upperHKSJ : (hk.ciHKSJ && hk.ciHKSJ[1]);
    if (isFinite(lo) && isFinite(hi)) {
      out.twoStageHKSJ = { est: re.pooled, lo, hi };
    }
  }
  return out;
}

function runExperiment({ k, nPerArm, alpha0, theta, sigmaAlpha, tau, nRep, seed0 }) {
  const acc = {
    oneStage: { est: [], cov: 0, n: 0, width: [] },
    twoStageZ: { est: [], cov: 0, n: 0, width: [] },
    twoStageHKSJ: { est: [], cov: 0, n: 0, width: [] },
  };
  for (let r = 0; r < nRep; r++) {
    const { rows } = generateIPD({ seed: seed0 + r * 7919, k, nPerArm, alpha0, theta, sigmaAlpha, tau });
    const f = fitOnce(rows);
    for (const key of ['oneStage', 'twoStageZ', 'twoStageHKSJ']) {
      const e = f[key];
      if (!e) continue;
      acc[key].n++;
      acc[key].est.push(e.est);
      acc[key].width.push(e.hi - e.lo);
      if (theta >= e.lo && theta <= e.hi) acc[key].cov++;
    }
  }
  const summ = {};
  for (const key of ['oneStage', 'twoStageZ', 'twoStageHKSJ']) {
    const a = acc[key];
    summ[key] = a.n === 0 ? null : {
      n: a.n, meanEst: mean(a.est), bias: mean(a.est) - theta,
      coverage: a.cov / a.n, meanWidth: mean(a.width),
    };
  }
  return { theta, tau, k, nPerArm, nRep, ...summ };
}

export { runExperiment, fitOnce, mean };

if (process.argv[1] && process.argv[1].endsWith('harness.mjs')) {
  const fmt = (s) => s ? `est=${s.meanEst.toFixed(3)} bias=${(s.bias>=0?'+':'')}${s.bias.toFixed(3)} cov=${(s.coverage*100).toFixed(1)}% w=${s.meanWidth.toFixed(2)} n=${s.n}` : 'n/a';
  console.log('IPD-META-PRO truth-recovery (binary, true log-OR)');
  console.log('one-stage = OneStageGLMM.fitBinary (PQL, RE study intercept)');
  console.log('two-stage = per-study logOR -> MetaAnalysis.randomEffectsREML (+HKSJ)\n');

  const scenarios = [
    { name: 'k=10 large, tau=0.1', k: 10, nPerArm: 300, alpha0: -1.0, theta: Math.log(0.7), sigmaAlpha: 0.5, tau: 0.1, nRep: 200 },
    { name: 'k=10 null,  tau=0.3', k: 10, nPerArm: 300, alpha0: -1.0, theta: 0,             sigmaAlpha: 0.5, tau: 0.3, nRep: 200 },
    { name: 'k=5  small, tau=0.3', k: 5,  nPerArm: 200, alpha0: -1.0, theta: Math.log(0.7), sigmaAlpha: 0.5, tau: 0.3, nRep: 200 },
  ];
  for (const sc of scenarios) {
    const res = runExperiment({ ...sc, seed0: 4242 });
    console.log(`-- ${sc.name}  (true log-OR=${sc.theta.toFixed(3)}, tau=${sc.tau}) --`);
    console.log(`   one-stage      : ${fmt(res.oneStage)}`);
    console.log(`   two-stage z    : ${fmt(res.twoStageZ)}`);
    console.log(`   two-stage HKSJ : ${fmt(res.twoStageHKSJ)}\n`);
  }
}
