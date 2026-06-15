/* node --test : truth-recovery validation for ipd-meta-pro.
 * Injects known-truth binary IPD (k studies, known overall log-OR, known
 * between-study heterogeneity tau, individual outcomes), runs the repo's OWN
 * one-stage GLMM and two-stage REML(+HKSJ) pooling, asserts bias + coverage. */
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { runExperiment, fitOnce } from './harness.mjs';
import { generateIPD, perStudyLogOR } from './dgp-ipd.mjs';
import { MetaAnalysis, OneStageGLMM } from './engine.mjs';

const big = { k: 10, nPerArm: 300, alpha0: -1.0, theta: Math.log(0.7), sigmaAlpha: 0.5, tau: 0.1, nRep: 150, seed0: 4242 };
const small = { k: 5, nPerArm: 200, alpha0: -1.0, theta: Math.log(0.7), sigmaAlpha: 0.5, tau: 0.3, nRep: 150, seed0: 4242 };

test('one-stage GLMM recovers the true log-OR (low bias)', () => {
  const r = runExperiment(big);
  assert.ok(r.oneStage, 'one-stage produced no result');
  assert.ok(Math.abs(r.oneStage.bias) < 0.05, `one-stage bias too large: ${r.oneStage.bias}`);
});

test('two-stage REML recovers the true log-OR (low bias)', () => {
  const r = runExperiment(big);
  assert.ok(r.twoStageZ, 'two-stage produced no result');
  assert.ok(Math.abs(r.twoStageZ.bias) < 0.05, `two-stage bias too large: ${r.twoStageZ.bias}`);
});

test('one-stage and two-stage point estimates agree closely (k=10)', () => {
  const r = runExperiment(big);
  const gap = Math.abs(r.oneStage.meanEst - r.twoStageZ.meanEst);
  assert.ok(gap < 0.03, `one- vs two-stage point estimates disagree: ${gap}`);
});

test('two-stage z achieves near-nominal coverage at k=10', () => {
  const r = runExperiment(big);
  assert.ok(r.twoStageZ.coverage >= 0.90 && r.twoStageZ.coverage <= 0.99,
    `two-stage z coverage off nominal: ${r.twoStageZ.coverage}`);
});

test('HKSJ repairs small-k under-coverage vs naive z', () => {
  // At k=5 with heterogeneity, naive z under-covers; HKSJ (t_{k-1}) restores it.
  const r = runExperiment(small);
  assert.ok(r.twoStageHKSJ.coverage > r.twoStageZ.coverage,
    `HKSJ did not widen coverage: HKSJ=${r.twoStageHKSJ.coverage} z=${r.twoStageZ.coverage}`);
  assert.ok(r.twoStageHKSJ.coverage >= 0.93,
    `HKSJ coverage still too low at k=5: ${r.twoStageHKSJ.coverage}`);
});

test('per-study 2x2 -> log-OR + Woolf variance are finite and sane', () => {
  const { rows } = generateIPD({ ...big });
  const { effects, variances } = perStudyLogOR(rows);
  assert.equal(effects.length, big.k);
  assert.ok(effects.every(isFinite) && variances.every(v => isFinite(v) && v > 0));
});

test('single replicate produces all three fits', () => {
  const { rows } = generateIPD({ ...big, seed: 999 });
  const f = fitOnce(rows);
  assert.ok(f.oneStage && f.twoStageZ && f.twoStageHKSJ, 'a fit path returned null');
});
