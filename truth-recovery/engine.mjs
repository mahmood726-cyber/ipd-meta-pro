/* ============================================================================
 * truth-recovery/engine.mjs
 * Pure statistical engine bundled VERBATIM from ipd-meta-pro dev modules:
 *   - dev/modules/02_06_stats.js          (Stats: distributions/quantiles)
 *   - dev/modules/02_07_confidence-utils.js (getConfZ)
 *   - dev/modules/02_08_meta-analysis.js  (MetaAnalysis: two-stage pooling)
 *   - dev/modules/02_25_one-stage-glmm.js (OneStageGLMM: one-stage PQL)
 * Additive only: source modules unchanged. Browser globals stubbed; exports
 * appended at the end. No source line edited.
 * ==========================================================================*/

// --- Minimal browser-global stubs (the modules feature-detect these) ---
const APP = { config: { confLevel: 0.95 } };
const window = undefined;  // modules guard with `typeof window !== 'undefined'`

// ===================== BEGIN 02_06_stats.js =====================
const Stats = {



 normalPDF: (x, mean = 0, sd = 1) => {

 const z = (x - mean) / sd;

 return Math.exp(-0.5 * z * z) / (sd * Math.sqrt(2 * Math.PI));

 },



 normalCDF: (x, mean = 0, sd = 1) => {

 const z = (x - mean) / sd;

 const t = 1 / (1 + 0.2316419 * Math.abs(z));

 const d = 0.3989422804 * Math.exp(-z * z / 2);

 const p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));

 return z > 0 ? 1 - p : p;

 },



 normalQuantile: (p, mean = 0, sd = 1) => {



 if (p <= 0) return -Infinity;

 if (p >= 1) return Infinity;

 if (p === 0.5) return mean;



 const a = [

 -3.969683028665376e+01, 2.209460984245205e+02,

 -2.759285104469687e+02, 1.383577518672690e+02,

 -3.066479806614716e+01, 2.506628277459239e+00

 ];

 const b = [

 -5.447609879822406e+01, 1.615858368580409e+02,

 -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01

 ];

 const c = [

 -7.784894002430293e-03, -3.223964580411365e-01,

 -2.400758277161838e+00, -2.549732539343734e+00,

 4.374664141464968e+00, 2.938163982698783e+00

 ];

 const d = [

 7.784695709041462e-03, 3.224671290700398e-01,

 2.445134137142996e+00, 3.754408661907416e+00

 ];



 const pLow = 0.02425;

 const pHigh = 1 - pLow;

 let q, r;



 if (p < pLow) {

 q = Math.sqrt(-2 * Math.log(p));

 return mean + sd * (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /

 ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);

 } else if (p <= pHigh) {

 q = p - 0.5;

 r = q * q;

 return mean + sd * (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /

 (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);

 } else {

 q = Math.sqrt(-2 * Math.log(1 - p));

 return mean + sd * -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /

 ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);

 }

 },



 chiSquareCDF: (x, df) => {

 if (x <= 0) return 0;

 return Stats.gammaCDF(x / 2, df / 2);

 },



 chiSquareQuantile: (p, df) => {



 if (p <= 0) return 0;

 if (p >= 1) return Infinity;



 let x = df;

 for (let i = 0; i < 100; i++) {

 const fx = Stats.chiSquareCDF(x, df) - p;

 const fpx = Stats.chiSquarePDF(x, df);

 if (Math.abs(fx) < 1e-10 || fpx === 0) break;

 x = Math.max(0.001, x - fx / fpx);

 }

 return x;

 },



 chiSquarePDF: (x, df) => {

 if (x <= 0) return 0;

 const k = df / 2;

 return Math.pow(x, k - 1) * Math.exp(-x / 2) / (Math.pow(2, k) * Stats.gamma(k));

 },



 gamma: (z) => {

 const g = 7;

 const c = [

 0.99999999999980993, 676.5203681218851, -1259.1392167224028,

 771.32342877765313, -176.61502916214059, 12.507343278686905,

 -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7

 ];



 if (z < 0.5) {

 return Math.PI / (Math.sin(Math.PI * z) * Stats.gamma(1 - z));

 }



 z -= 1;

 let x = c[0];

 for (let i = 1; i < g + 2; i++) {

 x += c[i] / (z + i);

 }

 const t = z + g + 0.5;

 return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;

 },



 gammaCDF: (x, a) => {



 if (x <= 0) return 0;

 if (x > 1000) return 1;



 const EPSILON = 1e-12;

 let sum = 0;

 let term = 1 / a;

 sum = term;



 for (let n = 1; n < 1000; n++) {

 term *= x / (a + n);

 sum += term;

 if (Math.abs(term) < EPSILON) break;

 }



 return Math.exp(-x + a * Math.log(x) - Math.log(Stats.gamma(a))) * sum;

 },



 tCDF: (t, df) => {

 const x = df / (df + t * t);

 return 1 - 0.5 * Stats.betaCDF(x, df / 2, 0.5);

 },



 tQuantile: (p, df) => {

 if (p === 0.5) return 0;

 const sign = p < 0.5 ? -1 : 1;

 const p2 = p < 0.5 ? p : 1 - p;



 const a = 1 / (df - 0.5);

 const b = 48 / (a * a);

 let c = ((20700 * a / b - 98) * a - 16) * a + 96.36;

 const d = ((94.5 / (b + c) - 3) / b + 1) * Math.sqrt(a * Math.PI / 2) * df;

 let x = d * p2;

 let y = Math.pow(x, 2 / df);



 if (y > 0.05 + a) {

 x = Stats.normalQuantile(p2);

 y = x * x;

 if (df < 5) c += 0.3 * (df - 4.5) * (x + 0.6);

 c += (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b;

 y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;

 y = a * y * y;

 y = y > 0.002 ? Math.exp(y) - 1 : 0.5 * y * y + y;

 } else {

 y = ((1 / (((df + 6) / (df * y) - 0.089 * d - 0.822) * (df + 2) * 3) + 0.5 / (df + 4)) * y - 1) * (df + 1) / (df + 2) + 1 / y;

 }



 return sign * Math.sqrt(df * y);

 },



 betaCDF: (x, a, b) => {

 if (x <= 0) return 0;

 if (x >= 1) return 1;



 const bt = Math.exp(

 Stats.logGamma(a + b) - Stats.logGamma(a) - Stats.logGamma(b) +

 a * Math.log(x) + b * Math.log(1 - x)

 );



 if (x < (a + 1) / (a + b + 2)) {

 return bt * Stats.betaCF(x, a, b) / a;

 } else {

 return 1 - bt * Stats.betaCF(1 - x, b, a) / b;

 }

 },



 betaCF: (x, a, b) => {

 const EPSILON = 1e-10;

 const MAX_ITER = 200;



 let qab = a + b;

 let qap = a + 1;

 let qam = a - 1;

 let c = 1;

 let d = 1 - qab * x / qap;

 if (Math.abs(d) < EPSILON) d = EPSILON;

 d = 1 / d;

 let h = d;



 for (let m = 1; m <= MAX_ITER; m++) {

 let m2 = 2 * m;

 let aa = m * (b - m) * x / ((qam + m2) * (a + m2));

 d = 1 + aa * d;

 if (Math.abs(d) < EPSILON) d = EPSILON;

 c = 1 + aa / c;

 if (Math.abs(c) < EPSILON) c = EPSILON;

 d = 1 / d;

 h *= d * c;



 aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

 d = 1 + aa * d;

 if (Math.abs(d) < EPSILON) d = EPSILON;

 c = 1 + aa / c;

 if (Math.abs(c) < EPSILON) c = EPSILON;

 d = 1 / d;

 const del = d * c;

 h *= del;



 if (Math.abs(del - 1) < EPSILON) break;

 }



 return h;

 },



 logGamma: (x) => {

 const c = [

 76.18009172947146, -86.50532032941677, 24.01409824083091,

 -1.231739572450155, 0.001208650973866179, -0.000005395239384953

 ];



 let y = x;

 let tmp = x + 5.5;

 tmp -= (x + 0.5) * Math.log(tmp);

 let ser = 1.000000000190015;



 for (let j = 0; j < 6; j++) {

 ser += c[j] / ++y;

 }



 return -tmp + Math.log(2.5066282746310005 * ser / x);

 },



 mean: (arr) => arr.reduce((a, b) => a + b, 0) / arr.length,



variance: (arr, ddof = 0) => {

if (!Array.isArray(arr) || arr.length === 0) return NaN;

const divisor = arr.length - ddof;

if (divisor <= 0) return NaN;

const m = Stats.mean(arr);

return arr.reduce((acc, val) => acc + Math.pow(val - m, 2), 0) / divisor;

},



sd: (arr, ddof = 0) => Math.sqrt(Stats.variance(arr, ddof)),



 median: (arr) => {

 const sorted = [...arr].sort((a, b) => a - b);

 const mid = Math.floor(sorted.length / 2);

 return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;

 },



 quantile: (arr, p) => {

 const sorted = [...arr].sort((a, b) => a - b);

 const pos = (sorted.length - 1) * p;

 const base = Math.floor(pos);

 const rest = pos - base;

 if (sorted[base + 1] !== undefined) {

 return sorted[base] + rest * (sorted[base + 1] - sorted[base]);

 }

 return sorted[base];

 },



 weightedMean: (values, weights) => {

 let sum = 0, weightSum = 0;

 for (let i = 0; i < values.length; i++) {

 sum += values[i] * weights[i];

 weightSum += weights[i];

 }

 return sum / weightSum;

 },



 linearRegression: (x, y, weights = null) => {

 const n = x.length;

 const w = weights || new Array(n).fill(1);



 let sumW = 0, sumWX = 0, sumWY = 0, sumWXX = 0, sumWXY = 0;



 for (let i = 0; i < n; i++) {

 sumW += w[i];

 sumWX += w[i] * x[i];

 sumWY += w[i] * y[i];

 sumWXX += w[i] * x[i] * x[i];

 sumWXY += w[i] * x[i] * y[i];

 }



 const slope = (sumW * sumWXY - sumWX * sumWY) / (sumW * sumWXX - sumWX * sumWX);

 const intercept = (sumWY - slope * sumWX) / sumW;



 let ssRes = 0, ssTot = 0;

 const yMean = sumWY / sumW;



 for (let i = 0; i < n; i++) {

 const yPred = intercept + slope * x[i];

 ssRes += w[i] * Math.pow(y[i] - yPred, 2);

 ssTot += w[i] * Math.pow(y[i] - yMean, 2);

 }



 const rSquared = 1 - ssRes / ssTot;

 const mse = ssRes / (sumW - 2);

 const seSlope = Math.sqrt(mse * sumW / (sumW * sumWXX - sumWX * sumWX));

 const seIntercept = Math.sqrt(mse * sumWXX / (sumW * sumWXX - sumWX * sumWX));



 return { slope, intercept, rSquared, seSlope, seIntercept, mse };

 }

 };



 // Confidence level z-value: uses APP.config.confLevel (default 0.95 → z=1.96)

 // Call getConfZ() instead of hardcoding 1.96 so CIs respect user's chosen level

 // ===================== END 02_06_stats.js =====================

// ===================== BEGIN 02_07_confidence-utils.js =====================
function getConfZ() {

 const level = (typeof APP !== 'undefined' && APP.config && APP.config.confLevel) ? APP.config.confLevel : 0.95;

 const alpha = 1 - level;

 // Use Stats.normalQuantile if available, otherwise approximate

 if (typeof Stats !== 'undefined' && Stats.normalQuantile) {

 return Stats.normalQuantile(1 - alpha / 2);

 }

 // Rational approximation fallback (Abramowitz & Stegun 26.2.23)

 const p = 1 - alpha / 2;

 const t = Math.sqrt(-2 * Math.log(1 - p));

 return t - (2.515517 + 0.802853 * t + 0.010328 * t * t) / (1 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t);

 }



 // ===================== END 02_07_confidence-utils.js =====================

// ===================== BEGIN 02_08_meta-analysis.js =====================
const MetaAnalysis = {

 MIN_VARIANCE: 1e-10,

 _prepareFE: (effects, variances) => {

 const weights = variances.map(v => 1 / Math.max(v, MetaAnalysis.MIN_VARIANCE));

 const sumW = weights.reduce((a, b) => a + b, 0);

 const pooled = effects.reduce((sum, e, i) => sum + weights[i] * e, 0) / sumW;

 let Q = 0;
 for (let i = 0; i < effects.length; i++) {
 Q += weights[i] * Math.pow(effects[i] - pooled, 2);
 }

 const sumW2 = weights.reduce((a, b) => a + b * b, 0);
 const C = sumW > 0 ? (sumW - sumW2 / sumW) : 0;
 const df = Math.max(0, effects.length - 1);
 const pQ = df > 0 ? (1 - Stats.chiSquareCDF(Q, df)) : 1;
 const I2Q = (df > 0 && Q > df && Q > 0) ? Math.max(0, 100 * (Q - df) / Q) : 0;
 const H2Q = df > 0 ? (Q / df) : 1;

 return { weights, sumW, sumW2, pooled, Q, C, df, pQ, I2Q, H2Q };

 },

 _withSummaryFields: (result, estimate = result.pooled, se = result.se, options = {}) => {

 const confLevel = (typeof APP !== 'undefined' && APP.config && Number.isFinite(APP.config.confLevel)) ? APP.config.confLevel : 0.95;
 const alpha = 1 - confLevel;
 const ciDf = Number.isFinite(options.ciDf) ? options.ciDf : null;
 const pValue = Number.isFinite(options.pValue) ? options.pValue :
  (Number.isFinite(result.pValue) ? result.pValue : (Number.isFinite(result.p) ? result.p : NaN));
 const k = Number.isFinite(options.k) ? options.k :
  (Number.isFinite(result.k) ? result.k : (Number.isFinite(result.df) ? result.df + 1 : undefined));

 let lower = NaN;
 let upper = NaN;

 if (Number.isFinite(estimate) && Number.isFinite(se)) {
  const crit = (ciDf !== null && ciDf > 0)
   ? Stats.tQuantile(1 - alpha / 2, ciDf)
   : getConfZ();
  lower = estimate - crit * se;
  upper = estimate + crit * se;
 }

 return {
  ...result,
  pooled: estimate,
  effect: estimate,
  estimate: estimate,
  lower: lower,
  upper: upper,
  ci: [lower, upper],
  ci_lower: lower,
  ci_upper: upper,
  pValue: pValue,
  k: k
 };

 },

 _buildRandomResult: (effects, variances, tau2) => {

 const tau2Safe = Math.max(0, Number(tau2) ?? 0);
 const fe = MetaAnalysis._prepareFE(effects, variances);

 const weights = variances.map(v => 1 / Math.max(v + tau2Safe, MetaAnalysis.MIN_VARIANCE));
 const sumW = weights.reduce((a, b) => a + b, 0);
 const pooled = effects.reduce((sum, e, i) => sum + weights[i] * e, 0) / sumW;
 const variance = 1 / sumW;
 const se = Math.sqrt(variance);
 const z = pooled / se;
 const p = 2 * (1 - Stats.normalCDF(Math.abs(z)));

 let I2 = 0;
 let H2 = 1;
 if (tau2Safe > 0 && fe.C > 0 && fe.df > 0) {
 const vTilde = fe.df / fe.C;
 if (vTilde > 0) {
 I2 = 100 * tau2Safe / (tau2Safe + vTilde);
 H2 = 1 / Math.max(MetaAnalysis.MIN_VARIANCE, 1 - I2 / 100);
 }
 }

 if (!isFinite(I2) || I2 < 0) I2 = 0;
 if (!isFinite(H2) || H2 < 1) H2 = 1;

 return MetaAnalysis._withSummaryFields({
 pooled, variance, se, z, p,
 weights,
 tau2: tau2Safe,
 tau: Math.sqrt(tau2Safe),
 Q: fe.Q,
 pQ: fe.pQ,
 I2: I2,
 H2: H2,
 df: fe.df
 }, pooled, se, { pValue: p, k: effects.length });

 },

 _estimateTau2PM: (effects, variances) => {

 const k = effects.length;
 const df = k - 1;
 if (k < 2) return 0;

 const qAt = (tau2) => {
 const weights = variances.map(v => 1 / Math.max(v + tau2, MetaAnalysis.MIN_VARIANCE));
 const sumW = weights.reduce((a, b) => a + b, 0);
 const pooled = effects.reduce((sum, e, i) => sum + weights[i] * e, 0) / sumW;
 let q = 0;
 for (let i = 0; i < k; i++) q += weights[i] * Math.pow(effects[i] - pooled, 2);
 return q;
 };

 if (qAt(0) <= df) return 0;

 let low = 0;
 let high = Math.max(1e-6, MetaAnalysis.randomEffectsDL(effects, variances).tau2 ?? 1e-6);
 while (qAt(high) > df && high < 1e6) high *= 2;

 for (let i = 0; i < 120; i++) {
 const mid = (low + high) / 2;
 const qMid = qAt(mid);
 if (qMid > df) low = mid;
 else high = mid;
 if (Math.abs(high - low) < 1e-10) break;
 }

 return Math.max(0, (low + high) / 2);

 },

 _weightedProjectionStats: (effects, variances, tau2) => {

 const weights = variances.map(v => 1 / Math.max(v + tau2, MetaAnalysis.MIN_VARIANCE));
 const sumW = weights.reduce((a, b) => a + b, 0);
 const sumW2 = weights.reduce((a, b) => a + b * b, 0);
 const sumW3 = weights.reduce((a, b) => a + b * b * b, 0);

 const pooled = effects.reduce((sum, e, i) => sum + weights[i] * e, 0) / sumW;

 let rss = 0;
 let yPPy = 0;
 for (let i = 0; i < effects.length; i++) {
 const resid = effects[i] - pooled;
 rss += weights[i] * resid * resid;
 const py = weights[i] * resid;
 yPPy += py * py;
 }

 const traceP = sumW - sumW2 / sumW;
 const tracePP = sumW2 - (2 * sumW3) / sumW + (sumW2 * sumW2) / (sumW * sumW);

 return { weights, sumW, sumW2, pooled, rss, yPPy, traceP, tracePP };

 },

 _estimateTau2Fisher: (effects, variances, useREML, options = {}) => {

 const k = effects.length;
 const p = 1;
 if (k <= p) return 0;

 const threshold = options.threshold ?? 1e-5;
 const maxIter = options.maxIter ?? 100;
 const tau2Min = options.tau2Min ?? 0;
 const tol = options.tol ?? 1e-4;
 const ll0check = options.ll0check ?? true;

 let tau2 = MetaAnalysis._estimateTau2HE(effects, variances);
 tau2 = Math.max(tau2Min, tau2);

 let change = threshold + 1;
 let iter = 0;

 while (change > threshold) {
 iter += 1;
 const oldTau2 = tau2;
 const stats = MetaAnalysis._weightedProjectionStats(effects, variances, tau2);

 let adj = 0;
 if (useREML) {
 if (stats.tracePP > 0) {
 adj = (stats.yPPy - stats.traceP) / stats.tracePP;
 }
 } else if (stats.sumW2 > 0) {
 adj = (stats.yPPy - stats.sumW) / stats.sumW2;
 }

 if (!isFinite(adj)) adj = 0;

 while (tau2 + adj < tau2Min) adj /= 2;
 tau2 += adj;
 change = Math.abs(oldTau2 - tau2);

 if (iter > maxIter) break;
 }

 tau2 = Math.max(tau2Min, tau2);

 if (ll0check && tau2 > threshold && variances.every(v => v > 0)) {
 const llStats0 = MetaAnalysis._weightedProjectionStats(effects, variances, 0);
 const llStats = MetaAnalysis._weightedProjectionStats(effects, variances, tau2);
 const kEff = k - p;

 let ll0;
 let ll;
 if (useREML) {
 ll0 = -0.5 * kEff * Math.log(2 * Math.PI) - 0.5 * variances.reduce((s, v) => s + Math.log(v), 0) - 0.5 * Math.log(Math.max(llStats0.sumW, MetaAnalysis.MIN_VARIANCE)) - 0.5 * llStats0.rss;
 ll = -0.5 * kEff * Math.log(2 * Math.PI) - 0.5 * variances.reduce((s, v) => s + Math.log(v + tau2), 0) - 0.5 * Math.log(Math.max(llStats.sumW, MetaAnalysis.MIN_VARIANCE)) - 0.5 * llStats.rss;
 } else {
 ll0 = -0.5 * k * Math.log(2 * Math.PI) - 0.5 * variances.reduce((s, v) => s + Math.log(v), 0) - 0.5 * llStats0.rss;
 ll = -0.5 * k * Math.log(2 * Math.PI) - 0.5 * variances.reduce((s, v) => s + Math.log(v + tau2), 0) - 0.5 * llStats.rss;
 }

 if (ll0 - ll > tol) tau2 = 0;
 }

 return Math.max(tau2Min, tau2);

 },

 _estimateTau2HE: (effects, variances) => {
 const k = effects.length;
 if (k < 2) return 0;
 const meanEffect = effects.reduce((a, b) => a + b, 0) / k;
 let ss = 0;
 for (let i = 0; i < k; i++) ss += Math.pow(effects[i] - meanEffect, 2);
 const sampleVar = ss / (k - 1);
 const meanVar = variances.reduce((a, b) => a + b, 0) / k;
 return Math.max(0, sampleVar - meanVar);
 },

 _estimateTau2SJ: (effects, variances) => {
 const k = effects.length;
 const p = 1;
 if (k <= p) return 0;

 const meanEffect = effects.reduce((a, b) => a + b, 0) / k;
 let ss = 0;
 for (let i = 0; i < k; i++) ss += Math.pow(effects[i] - meanEffect, 2);
 const tau2_0 = ss / k;

 const stats = MetaAnalysis._weightedProjectionStats(effects, variances, tau2_0);
 return Math.max(0, tau2_0 * stats.rss / (k - p));
 },

 fixedEffect: (effects, variances) => {

 const fe = MetaAnalysis._prepareFE(effects, variances);
 const variance = 1 / fe.sumW;
 const se = Math.sqrt(variance);
 const z = fe.pooled / se;
 const p = 2 * (1 - Stats.normalCDF(Math.abs(z)));

 return MetaAnalysis._withSummaryFields({
 pooled: fe.pooled,
 variance: variance,
 se: se,
 z: z,
 p: p,
 weights: fe.weights,
 tau2: 0,
 tau: 0,
 Q: fe.Q,
 pQ: fe.pQ,
 I2: fe.I2Q,
 H2: fe.H2Q,
 df: fe.df
 }, fe.pooled, se, { pValue: p, k: effects.length });

 },

 randomEffectsDL: (effects, variances) => {
 const fe = MetaAnalysis._prepareFE(effects, variances);
 const tau2 = fe.C > 0 ? Math.max(0, (fe.Q - fe.df) / fe.C) : 0;
 return MetaAnalysis._buildRandomResult(effects, variances, tau2);
 },

 randomEffectsREML: (effects, variances) => {
 const tau2 = MetaAnalysis._estimateTau2Fisher(effects, variances, true);
 return MetaAnalysis._buildRandomResult(effects, variances, tau2);
 },

 randomEffectsPM: (effects, variances) => {
 const tau2 = MetaAnalysis._estimateTau2PM(effects, variances);
 return MetaAnalysis._buildRandomResult(effects, variances, tau2);
 },

 randomEffectsSJ: (effects, variances) => {
 const tau2 = MetaAnalysis._estimateTau2SJ(effects, variances);
 return MetaAnalysis._buildRandomResult(effects, variances, tau2);
 },

 randomEffectsHE: (effects, variances) => {
 const tau2 = MetaAnalysis._estimateTau2HE(effects, variances);
 return MetaAnalysis._buildRandomResult(effects, variances, tau2);
 },

 randomEffectsML: (effects, variances) => {
 const tau2 = MetaAnalysis._estimateTau2Fisher(effects, variances, false);
 return MetaAnalysis._buildRandomResult(effects, variances, tau2);
 },



 applyHKSJ: (result, effects, variances) => {

 const k = effects.length;

 const weights = result.weights;

 const sumW = weights.reduce((a, b) => a + b, 0);



 let qAdj = 0;

 for (let i = 0; i < k; i++) {

 qAdj += weights[i] * Math.pow(effects[i] - result.pooled, 2);

 }

 qAdj /= (k - 1);



 const seAdj = result.se * Math.sqrt(Math.max(1, qAdj));

 const t = result.pooled / seAdj;

 const p = 2 * (1 - Stats.tCDF(Math.abs(t), k - 1));



 const adjusted = MetaAnalysis._withSummaryFields({
 ...result,
 seHKSJ: seAdj,
 tHKSJ: t,
 pHKSJ: p
 }, result.pooled, seAdj, { ciDf: k - 1, pValue: p, k: k });

 return {
 ...adjusted,
 seUnadjusted: result.se,
 pUnadjusted: result.p,
 lowerUnadjusted: result.lower,
 upperUnadjusted: result.upper,
 ciUnadjusted: result.ci,
 lowerHKSJ: adjusted.lower,
 upperHKSJ: adjusted.upper,
 ciHKSJ: adjusted.ci
 };

 },



 predictionInterval: (result, confLevel = 0.95) => {

 const k = result.df + 1;

 const tCrit = Stats.tQuantile(1 - (1 - confLevel) / 2, k - 2);

 const piWidth = tCrit * Math.sqrt(result.variance + result.tau2);



 return {

 lower: result.pooled - piWidth,

 upper: result.pooled + piWidth

 };

 },



 confidenceInterval: (estimate, se, confLevel = 0.95, df = null) => {

 const alpha = 1 - confLevel;

 let crit;

 if (df) {

 crit = Stats.tQuantile(1 - alpha / 2, df);

 } else {

 crit = Stats.normalQuantile(1 - alpha / 2);

 }

 return {

 lower: estimate - crit * se,

 upper: estimate + crit * se

 };

 }

 };



 // ===================== END 02_08_meta-analysis.js =====================

// ===================== BEGIN 02_25_one-stage-glmm.js =====================
// ═══════════════════════════════════════════════════════════════════
// One-Stage GLMM Engine — Binary outcomes via PQL (Penalized Quasi-Likelihood)
// Reference: Breslow & Clayton (1993), Debray et al. (2015, 2017)
// ═══════════════════════════════════════════════════════════════════

var OneStageGLMM = (function() {
    'use strict';

    // --- Matrix helpers (small-scale, for p×p where p is typically 2-10) ---

    function matCreate(rows, cols, fill) {
        var m = new Array(rows);
        for (var i = 0; i < rows; i++) {
            m[i] = new Array(cols);
            for (var j = 0; j < cols; j++) m[i][j] = (fill !== undefined ? fill : 0);
        }
        return m;
    }

    function matClone(A) {
        return A.map(function(row) { return row.slice(); });
    }

    function matTranspose(A) {
        var m = A.length, n = A[0].length;
        var T = matCreate(n, m);
        for (var i = 0; i < m; i++)
            for (var j = 0; j < n; j++) T[j][i] = A[i][j];
        return T;
    }

    function matMul(A, B) {
        var m = A.length, n = B[0].length, k = B.length;
        var C = matCreate(m, n);
        for (var i = 0; i < m; i++)
            for (var j = 0; j < n; j++) {
                var s = 0;
                for (var l = 0; l < k; l++) s += A[i][l] * B[l][j];
                C[i][j] = s;
            }
        return C;
    }

    function matVecMul(A, v) {
        var n = A.length;
        var result = new Array(n);
        for (var i = 0; i < n; i++) {
            var s = 0;
            for (var j = 0; j < A[i].length; j++) s += A[i][j] * v[j];
            result[i] = s;
        }
        return result;
    }

    function matAdd(A, B) {
        var m = A.length, n = A[0].length;
        var C = matCreate(m, n);
        for (var i = 0; i < m; i++)
            for (var j = 0; j < n; j++) C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    // Cholesky decomposition: A = L L^T, returns L (lower triangular)
    function cholesky(A) {
        var n = A.length;
        var L = matCreate(n, n);
        for (var i = 0; i < n; i++) {
            for (var j = 0; j <= i; j++) {
                var s = 0;
                for (var k = 0; k < j; k++) s += L[i][k] * L[j][k];
                if (i === j) {
                    var val = A[i][i] - s;
                    if (val <= 0) val = 1e-12; // Regularize
                    L[i][j] = Math.sqrt(val);
                } else {
                    L[i][j] = (A[i][j] - s) / L[j][j];
                }
            }
        }
        return L;
    }

    // Solve L x = b (forward substitution) where L is lower triangular
    function forwardSolve(L, b) {
        var n = L.length;
        var x = new Array(n);
        for (var i = 0; i < n; i++) {
            var s = b[i];
            for (var j = 0; j < i; j++) s -= L[i][j] * x[j];
            x[i] = s / L[i][i];
        }
        return x;
    }

    // Solve L^T x = b (backward substitution) where L is lower triangular
    function backwardSolve(L, b) {
        var n = L.length;
        var x = new Array(n);
        for (var i = n - 1; i >= 0; i--) {
            var s = b[i];
            for (var j = i + 1; j < n; j++) s -= L[j][i] * x[j];
            x[i] = s / L[i][i];
        }
        return x;
    }

    // Solve A x = b using Cholesky: A = L L^T, then L y = b, L^T x = y
    function cholSolve(A, b) {
        var L = cholesky(A);
        var y = forwardSolve(L, b);
        return backwardSolve(L, y);
    }

    // Invert symmetric positive definite matrix via Cholesky
    function cholInverse(A) {
        var n = A.length;
        var inv = matCreate(n, n);
        var L = cholesky(A);
        for (var col = 0; col < n; col++) {
            var e = new Array(n).fill(0);
            e[col] = 1;
            var y = forwardSolve(L, e);
            var x = backwardSolve(L, y);
            for (var row = 0; row < n; row++) inv[row][col] = x[row];
        }
        return inv;
    }

    // Log determinant via Cholesky: log|A| = 2 * sum(log(L_ii))
    function cholLogDet(A) {
        var L = cholesky(A);
        var s = 0;
        for (var i = 0; i < L.length; i++) s += Math.log(Math.max(L[i][i], 1e-300));
        return 2 * s;
    }

    // --- Link functions ---

    var LINKS = {
        logit: {
            link: function(mu) { return Math.log(mu / (1 - mu)); },
            invLink: function(eta) {
                eta = Math.max(-500, Math.min(500, eta));
                return 1 / (1 + Math.exp(-eta));
            },
            dmu: function(mu) { return mu * (1 - mu); }, // d(mu)/d(eta) = variance function for logistic
            variance: function(mu) { return mu * (1 - mu); }
        },
        log: {
            link: function(mu) { return Math.log(Math.max(mu, 1e-300)); },
            invLink: function(eta) { return Math.exp(Math.min(eta, 500)); },
            dmu: function(mu) { return mu; },
            variance: function(mu) { return mu; }
        },
        identity: {
            link: function(mu) { return mu; },
            invLink: function(eta) { return eta; },
            dmu: function() { return 1; },
            variance: function() { return 1; }
        }
    };

    // --- Core PQL Algorithm ---
    // Fits: g(E[Y_ij | b_i]) = X_ij' β + Z_ij' b_i
    //       b_i ~ N(0, D)
    // Via iterative weighted LMM

    function fitPQL(y, X, studyIdx, options) {
        options = options || {};
        var family = options.family || 'binomial';
        var linkName = options.link || (family === 'binomial' ? 'logit' : 'identity');
        var linkFn = LINKS[linkName] || LINKS.logit;
        var maxOuterIter = options.maxIter || 50;
        var outerTol = options.tolerance || 1e-6;
        var verbose = options.verbose || false;

        var N = y.length;
        var p = X[0].length;

        // Identify unique studies
        var studyLabels = [];
        var studyMap = {};
        for (var i = 0; i < N; i++) {
            var sid = studyIdx[i];
            if (!studyMap.hasOwnProperty(sid)) {
                studyMap[sid] = studyLabels.length;
                studyLabels.push(sid);
            }
        }
        var K = studyLabels.length;

        // Study membership: studyOf[i] = index into studyLabels
        var studyOf = new Array(N);
        for (var ii = 0; ii < N; ii++) studyOf[ii] = studyMap[studyIdx[ii]];

        // Initialize: β from fixed-effect GLM (IRLS), b = 0
        var beta = new Array(p).fill(0);
        var b = new Array(K).fill(0); // Random intercepts
        var tau2 = 0.1; // Between-study variance

        // Initial IRLS for β (ignore random effects)
        for (var irls = 0; irls < 25; irls++) {
            var W_diag = new Array(N);
            var z_work = new Array(N);
            for (var j = 0; j < N; j++) {
                var eta_j = 0;
                for (var l = 0; l < p; l++) eta_j += X[j][l] * beta[l];
                var mu_j = linkFn.invLink(eta_j);
                if (family === 'binomial') mu_j = Math.max(1e-6, Math.min(1 - 1e-6, mu_j));
                var dmu = linkFn.dmu(mu_j);
                var v = linkFn.variance(mu_j);
                W_diag[j] = dmu * dmu / Math.max(v, 1e-10);
                z_work[j] = eta_j + (y[j] - mu_j) / Math.max(dmu, 1e-10);
            }
            // Weighted normal equations: (X'WX) β = X'Wz
            var XtWX = matCreate(p, p);
            var XtWz = new Array(p).fill(0);
            for (var r = 0; r < N; r++) {
                for (var a = 0; a < p; a++) {
                    XtWz[a] += X[r][a] * W_diag[r] * z_work[r];
                    for (var c = 0; c < p; c++) XtWX[a][c] += X[r][a] * W_diag[r] * X[r][c];
                }
            }
            // Regularize diagonal
            for (var d = 0; d < p; d++) XtWX[d][d] += 1e-8;
            var newBeta;
            try { newBeta = cholSolve(XtWX, XtWz); } catch (e) { break; }
            var maxDelta = 0;
            for (var dd = 0; dd < p; dd++) maxDelta = Math.max(maxDelta, Math.abs(newBeta[dd] - beta[dd]));
            beta = newBeta;
            if (maxDelta < 1e-8) break;
        }

        // --- PQL outer loop ---
        var converged = false;
        var outerIter = 0;
        var logLik = -Infinity;

        for (outerIter = 0; outerIter < maxOuterIter; outerIter++) {
            // E-step: compute working response and weights
            var W = new Array(N);
            var zw = new Array(N);
            for (var i2 = 0; i2 < N; i2++) {
                var eta = 0;
                for (var l2 = 0; l2 < p; l2++) eta += X[i2][l2] * beta[l2];
                eta += b[studyOf[i2]]; // Add random intercept
                var mu = linkFn.invLink(eta);
                if (family === 'binomial') mu = Math.max(1e-6, Math.min(1 - 1e-6, mu));
                var dm = linkFn.dmu(mu);
                var va = linkFn.variance(mu);
                W[i2] = dm * dm / Math.max(va, 1e-10);
                zw[i2] = eta + (y[i2] - mu) / Math.max(dm, 1e-10);
            }

            // M-step: fit weighted LMM on working response
            // (X'WX + Z'WZ/tau2)^{-1} [X'Wz; Z'Wz/tau2... no]
            // Henderson's mixed model equations:
            // [ X'WX      X'WZ    ] [β] = [X'Wz  ]
            // [ Z'WX   Z'WZ + D⁻¹ ] [b]   [Z'Wz  ]
            // where D = tau2 * I_K

            var pTotal = p + K;
            var LHS = matCreate(pTotal, pTotal);
            var RHS = new Array(pTotal).fill(0);

            // Build X'WX block and X'Wz
            for (var i3 = 0; i3 < N; i3++) {
                var wi = W[i3];
                var si = studyOf[i3];
                // X'WX
                for (var a2 = 0; a2 < p; a2++) {
                    RHS[a2] += X[i3][a2] * wi * zw[i3];
                    for (var c2 = 0; c2 < p; c2++) LHS[a2][c2] += X[i3][a2] * wi * X[i3][c2];
                    // X'WZ
                    LHS[a2][p + si] += X[i3][a2] * wi;
                    LHS[p + si][a2] += wi * X[i3][a2];
                }
                // Z'WZ
                LHS[p + si][p + si] += wi;
                // Z'Wz
                RHS[p + si] += wi * zw[i3];
            }

            // Add D^{-1} = (1/tau2) I_K to Z'WZ block
            var tau2Inv = 1 / Math.max(tau2, 1e-10);
            for (var k = 0; k < K; k++) {
                LHS[p + k][p + k] += tau2Inv;
            }

            // Regularize
            for (var d2 = 0; d2 < pTotal; d2++) LHS[d2][d2] += 1e-10;

            // Solve Henderson's equations
            var solution;
            try { solution = cholSolve(LHS, RHS); } catch (e) {
                // Fallback: increase regularization
                for (var d3 = 0; d3 < pTotal; d3++) LHS[d3][d3] += 1e-4;
                try { solution = cholSolve(LHS, RHS); } catch (e2) { break; }
            }

            var newBeta2 = solution.slice(0, p);
            var newB = solution.slice(p, p + K);

            // Update tau2 (REML-style)
            // tau2 = (b'b + trace(C22)) / K where C22 is posterior covariance of b
            var bSS = 0;
            for (var k2 = 0; k2 < K; k2++) bSS += newB[k2] * newB[k2];

            // Approximate trace of posterior covariance from diagonal of inverse
            // C22 = (Z'WZ + D^{-1})^{-1}
            var traceC22 = 0;
            try {
                var C22 = matCreate(K, K);
                for (var i4 = 0; i4 < K; i4++)
                    for (var j4 = 0; j4 < K; j4++)
                        C22[i4][j4] = LHS[p + i4][p + j4];
                var C22inv = cholInverse(C22);
                for (var k3 = 0; k3 < K; k3++) traceC22 += C22inv[k3][k3];
            } catch (e) {
                traceC22 = K * tau2 * 0.5; // Fallback
            }

            var newTau2 = Math.max(1e-8, (bSS + traceC22) / K);

            // Check convergence
            var maxD = 0;
            for (var dd2 = 0; dd2 < p; dd2++) maxD = Math.max(maxD, Math.abs(newBeta2[dd2] - beta[dd2]));
            for (var kk = 0; kk < K; kk++) maxD = Math.max(maxD, Math.abs(newB[kk] - b[kk]));
            maxD = Math.max(maxD, Math.abs(newTau2 - tau2));

            beta = newBeta2;
            b = newB;
            tau2 = newTau2;

            if (maxD < outerTol) {
                converged = true;
                break;
            }
        }

        // --- Compute standard errors from (X'WX)^{-1} ---
        // Recompute final W
        var Wfinal = new Array(N);
        for (var i5 = 0; i5 < N; i5++) {
            var eta5 = 0;
            for (var l5 = 0; l5 < p; l5++) eta5 += X[i5][l5] * beta[l5];
            eta5 += b[studyOf[i5]];
            var mu5 = linkFn.invLink(eta5);
            if (family === 'binomial') mu5 = Math.max(1e-6, Math.min(1 - 1e-6, mu5));
            Wfinal[i5] = linkFn.dmu(mu5) * linkFn.dmu(mu5) / Math.max(linkFn.variance(mu5), 1e-10);
        }

        // Marginal covariance of β: (X' V^{-1} X)^{-1} where V = W^{-1} + Z D Z'
        // For random intercepts: V_i = diag(1/w_ij) + tau2 * J
        // Using Woodbury: V_i^{-1} = W_i - W_i Z_i (Z_i' W_i Z_i + D^{-1})^{-1} Z_i' W_i
        var XtVinvX = matCreate(p, p);
        for (var s = 0; s < K; s++) {
            // Collect patients in this study
            var idx = [];
            for (var i6 = 0; i6 < N; i6++) if (studyOf[i6] === s) idx.push(i6);
            var nk = idx.length;

            // Z_k' W_k Z_k + D^{-1} = sum(w_ij) + 1/tau2
            var sumW = 0;
            for (var i7 = 0; i7 < nk; i7++) sumW += Wfinal[idx[i7]];
            var c_inv = sumW + tau2Inv;
            var c_val = 1 / c_inv;

            // X_k' V_k^{-1} X_k = X_k' W_k X_k - (X_k' W_k 1)(1' W_k X_k) * c
            for (var a3 = 0; a3 < p; a3++) {
                var xwSum_a = 0;
                for (var i8 = 0; i8 < nk; i8++) xwSum_a += X[idx[i8]][a3] * Wfinal[idx[i8]];
                for (var c3 = 0; c3 < p; c3++) {
                    var xwx = 0;
                    var xwSum_c = 0;
                    for (var i9 = 0; i9 < nk; i9++) {
                        xwx += X[idx[i9]][a3] * Wfinal[idx[i9]] * X[idx[i9]][c3];
                        xwSum_c += X[idx[i9]][c3] * Wfinal[idx[i9]];
                    }
                    XtVinvX[a3][c3] += xwx - xwSum_a * xwSum_c * c_val;
                }
            }
        }

        // Regularize and invert
        for (var d4 = 0; d4 < p; d4++) XtVinvX[d4][d4] += 1e-10;
        var betaCov;
        try { betaCov = cholInverse(XtVinvX); } catch (e) { betaCov = null; }
        var betaSE = new Array(p);
        if (betaCov) {
            for (var d5 = 0; d5 < p; d5++) betaSE[d5] = Math.sqrt(Math.max(0, betaCov[d5][d5]));
        } else {
            for (var d6 = 0; d6 < p; d6++) betaSE[d6] = NaN;
        }

        // Approximate log-likelihood (Laplace approximation to marginal)
        logLik = 0;
        for (var i10 = 0; i10 < N; i10++) {
            var eta10 = 0;
            for (var l10 = 0; l10 < p; l10++) eta10 += X[i10][l10] * beta[l10];
            eta10 += b[studyOf[i10]];
            var mu10 = linkFn.invLink(eta10);
            if (family === 'binomial') {
                mu10 = Math.max(1e-300, Math.min(1 - 1e-300, mu10));
                logLik += y[i10] * Math.log(mu10) + (1 - y[i10]) * Math.log(1 - mu10);
            } else {
                // Gaussian log-likelihood
                var resid10 = y[i10] - mu10;
                logLik += -0.5 * resid10 * resid10; // Up to constant
            }
        }
        // Penalty from random effects
        for (var k4 = 0; k4 < K; k4++) {
            logLik -= 0.5 * b[k4] * b[k4] / Math.max(tau2, 1e-10);
            logLik -= 0.5 * Math.log(2 * Math.PI * Math.max(tau2, 1e-10));
        }

        return {
            beta: beta,
            se: betaSE,
            vcov: betaCov,
            randomEffects: b,
            tau2: tau2,
            tau: Math.sqrt(tau2),
            sigma2: null, // Not applicable for binary
            nStudies: K,
            nPatients: N,
            nFixed: p,
            studyLabels: studyLabels,
            convergence: { converged: converged, iterations: outerIter + 1 },
            logLik: logLik,
            family: family,
            link: linkName
        };
    }

    // --- Public API ---

    return {

        /**
         * Fit a one-stage GLMM for binary outcomes (logistic mixed model)
         * @param {Array} data - IPD array of objects
         * @param {string} outcomeVar - Binary outcome variable name (0/1)
         * @param {string} treatmentVar - Treatment variable name (0/1)
         * @param {string} studyVar - Study identifier variable name
         * @param {Array} covariates - Optional covariate variable names
         * @param {Object} options - {family, link, maxIter, tolerance, interaction}
         * @returns {Object} Model fit results
         */
        fitBinary: function(data, outcomeVar, treatmentVar, studyVar, covariates, options) {
            covariates = covariates || [];
            options = options || {};

            // Clean and structure data
            var cleaned = [];
            (data || []).forEach(function(d) {
                var yRaw = d[outcomeVar];
                var y = (yRaw === 1 || yRaw === '1' || yRaw === true) ? 1 : (yRaw === 0 || yRaw === '0' || yRaw === false) ? 0 : null;
                if (y === null) return;
                var t = (d[treatmentVar] === 1 || d[treatmentVar] === '1' || d[treatmentVar] === true) ? 1 : 0;
                var sid = d[studyVar];
                if (sid === null || sid === undefined || sid === '') return;

                // Build covariate values
                var covVals = [];
                for (var c = 0; c < covariates.length; c++) {
                    var v = Number(d[covariates[c]]);
                    if (!isFinite(v)) return;
                    covVals.push(v);
                }

                // Design matrix row: [1, treatment, cov1, cov2, ..., interaction?]
                var x = [1, t];
                covVals.forEach(function(v) { x.push(v); });
                if (options.interaction && covariates.length > 0) {
                    // Add treatment × covariate interactions
                    covVals.forEach(function(v) { x.push(t * v); });
                }

                cleaned.push({ y: y, x: x, study: String(sid) });
            });

            if (cleaned.length < 20) {
                return { error: 'Need at least 20 complete cases for GLMM. Got ' + cleaned.length };
            }

            // Extract arrays
            var yArr = cleaned.map(function(r) { return r.y; });
            var XArr = cleaned.map(function(r) { return r.x; });
            var studyArr = cleaned.map(function(r) { return r.study; });

            // Fit PQL
            var fit = fitPQL(yArr, XArr, studyArr, {
                family: 'binomial',
                link: 'logit',
                maxIter: options.maxIter || 50,
                tolerance: options.tolerance || 1e-6
            });

            if (!fit || !fit.beta) return { error: 'GLMM failed to converge' };

            // Format results
            var treatIdx = 1;
            var treatEffect = fit.beta[treatIdx]; // log-OR
            var treatSE = fit.se[treatIdx];
            var z = typeof getConfZ === 'function' ? getConfZ() : 1.96;
            var OR = Math.exp(treatEffect);
            var OR_lower = Math.exp(treatEffect - z * treatSE);
            var OR_upper = Math.exp(treatEffect + z * treatSE);
            var zStat = treatEffect / treatSE;
            var pValue = 2 * (1 - (typeof Stats !== 'undefined' ? Stats.normalCDF(Math.abs(zStat)) : 0.5));

            // I² for treatment effect heterogeneity
            var tau2 = fit.tau2;
            var typicalSE2 = treatSE * treatSE;
            var I2 = (tau2 + typicalSE2) > 0 ? Math.max(0, Math.min(100, 100 * tau2 / (tau2 + typicalSE2))) : 0;

            // Prediction interval on log-OR scale
            var predLower = treatEffect - z * Math.sqrt(treatSE * treatSE + tau2);
            var predUpper = treatEffect + z * Math.sqrt(treatSE * treatSE + tau2);

            var result = {
                method: 'One-Stage GLMM (PQL)',
                family: 'binomial',
                link: 'logit',
                estimand: 'log-OR',
                nStudies: fit.nStudies,
                nPatients: fit.nPatients,
                studyLabels: fit.studyLabels,

                fixedEffects: {
                    names: ['intercept', 'treatment'].concat(covariates).concat(
                        options.interaction ? covariates.map(function(c) { return 'treatment:' + c; }) : []
                    ),
                    beta: fit.beta,
                    se: fit.se,
                    vcov: fit.vcov
                },

                treatment: {
                    logOR: treatEffect,
                    OR: OR,
                    se: treatSE,
                    CI: [OR_lower, OR_upper],
                    CI_logScale: [treatEffect - z * treatSE, treatEffect + z * treatSE],
                    predictionInterval: [Math.exp(predLower), Math.exp(predUpper)],
                    predictionInterval_logScale: [predLower, predUpper],
                    zStat: zStat,
                    pValue: pValue
                },

                randomEffects: {
                    tau2: tau2,
                    tau: fit.tau,
                    studyIntercepts: fit.randomEffects,
                    I2: I2
                },

                heterogeneity: {
                    tau_sq: tau2,
                    tau: fit.tau,
                    I2: I2
                },

                convergence: fit.convergence,
                logLik: fit.logLik,

                // Compatibility fields for two-stage pipeline
                pooled_effect: treatEffect,
                SE: treatSE,
                CI_lower: treatEffect - z * treatSE,
                CI_upper: treatEffect + z * treatSE,
                p_value: pValue,
                tau2: tau2,
                I2: I2
            };

            // Interaction effects (if covariates and interaction requested)
            if (options.interaction && covariates.length > 0) {
                var interactionResults = [];
                var baseP = 2 + covariates.length;
                for (var ci = 0; ci < covariates.length; ci++) {
                    var intIdx = baseP + ci;
                    if (intIdx < fit.beta.length) {
                        var intBeta = fit.beta[intIdx];
                        var intSE = fit.se[intIdx];
                        var intZ = intBeta / intSE;
                        var intP = 2 * (1 - (typeof Stats !== 'undefined' ? Stats.normalCDF(Math.abs(intZ)) : 0.5));
                        interactionResults.push({
                            covariate: covariates[ci],
                            beta: intBeta,
                            OR: Math.exp(intBeta),
                            se: intSE,
                            z: intZ,
                            p: intP
                        });
                    }
                }
                result.interactions = interactionResults;
            }

            return result;
        },

        /**
         * Fit a one-stage LMM for continuous outcomes with REML
         * Enhances existing runOneStageIPDMA with proper REML adjustment
         * @param {Array} data - IPD array of objects
         * @param {string} outcomeVar - Continuous outcome variable name
         * @param {string} treatmentVar - Treatment variable name (0/1)
         * @param {string} studyVar - Study identifier variable name
         * @param {Array} covariates - Optional covariate variable names
         * @param {Object} options - {maxIter, tolerance, randomSlope, interaction}
         * @returns {Object} Model fit results
         */
        fitContinuous: function(data, outcomeVar, treatmentVar, studyVar, covariates, options) {
            options = options || {};

            // If the existing runOneStageIPDMA is available, use it (it's already comprehensive)
            if (typeof runOneStageIPDMA === 'function' && !options.forceNew) {
                var existingResult = runOneStageIPDMA(data, outcomeVar, treatmentVar, studyVar, covariates);
                if (existingResult && !existingResult.error) {
                    existingResult.method = 'One-Stage LMM (Random Intercept + Slope)';
                    return existingResult;
                }
            }

            // Fallback: fit via PQL with identity link (equivalent to LMM)
            var cleaned = [];
            (data || []).forEach(function(d) {
                var y = Number(d[outcomeVar]);
                if (!isFinite(y)) return;
                var t = (d[treatmentVar] === 1 || d[treatmentVar] === '1') ? 1 : 0;
                var sid = d[studyVar];
                if (sid === null || sid === undefined || sid === '') return;

                var covVals = [];
                for (var c = 0; c < (covariates || []).length; c++) {
                    var v = Number(d[covariates[c]]);
                    if (!isFinite(v)) return;
                    covVals.push(v);
                }

                var x = [1, t];
                covVals.forEach(function(v) { x.push(v); });
                if (options.interaction && covariates && covariates.length > 0) {
                    covVals.forEach(function(v) { x.push(t * v); });
                }

                cleaned.push({ y: y, x: x, study: String(sid) });
            });

            if (cleaned.length < 10) {
                return { error: 'Need at least 10 complete cases. Got ' + cleaned.length };
            }

            var yArr = cleaned.map(function(r) { return r.y; });
            var XArr = cleaned.map(function(r) { return r.x; });
            var studyArr = cleaned.map(function(r) { return r.study; });

            var fit = fitPQL(yArr, XArr, studyArr, {
                family: 'gaussian',
                link: 'identity',
                maxIter: options.maxIter || 100,
                tolerance: options.tolerance || 1e-7
            });

            if (!fit || !fit.beta) return { error: 'LMM failed to converge' };

            var treatEffect = fit.beta[1];
            var treatSE = fit.se[1];
            var z = typeof getConfZ === 'function' ? getConfZ() : 1.96;

            return {
                method: 'One-Stage LMM (REML via PQL)',
                nStudies: fit.nStudies,
                nPatients: fit.nPatients,
                studyLabels: fit.studyLabels,
                fixedEffects: {
                    intercept: fit.beta[0],
                    treatment: treatEffect,
                    covariates: fit.beta.slice(2),
                    se: fit.se
                },
                treatment: {
                    effect: treatEffect,
                    se: treatSE,
                    CI: [treatEffect - z * treatSE, treatEffect + z * treatSE],
                    predictionInterval: [
                        treatEffect - z * Math.sqrt(treatSE * treatSE + fit.tau2),
                        treatEffect + z * Math.sqrt(treatSE * treatSE + fit.tau2)
                    ],
                    zStat: treatEffect / treatSE,
                    pValue: 2 * (1 - (typeof Stats !== 'undefined' ? Stats.normalCDF(Math.abs(treatEffect / treatSE)) : 0.5))
                },
                randomEffects: {
                    tau2: fit.tau2,
                    tau: fit.tau,
                    studyIntercepts: fit.randomEffects
                },
                heterogeneity: {
                    tau_sq: fit.tau2,
                    tau: fit.tau,
                    I2: fit.tau2 > 0 ? Math.max(0, Math.min(100, 100 * fit.tau2 / (fit.tau2 + treatSE * treatSE))) : 0
                },
                convergence: fit.convergence,
                logLik: fit.logLik,
                // Compatibility
                pooled_effect: treatEffect,
                SE: treatSE,
                CI_lower: treatEffect - z * treatSE,
                CI_upper: treatEffect + z * treatSE,
                p_value: 2 * (1 - (typeof Stats !== 'undefined' ? Stats.normalCDF(Math.abs(treatEffect / treatSE)) : 0.5)),
                tau2: fit.tau2,
                I2: fit.tau2 > 0 ? Math.max(0, Math.min(100, 100 * fit.tau2 / (fit.tau2 + treatSE * treatSE))) : 0
            };
        },

        // Expose internals for testing
        _fitPQL: fitPQL,
        _cholesky: cholesky,
        _cholSolve: cholSolve,
        _cholInverse: cholInverse,
        _matMul: matMul
    };
})();

if (typeof window !== 'undefined') {
    window.OneStageGLMM = OneStageGLMM;
}
// ===================== END 02_25_one-stage-glmm.js =====================

export { Stats, MetaAnalysis, OneStageGLMM, getConfZ };
