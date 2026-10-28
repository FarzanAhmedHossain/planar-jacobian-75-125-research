// Exact rational replay of the first two tau-faces of the distance-36 source PDE.
// This checks a compatible local jet only; it is not a full Keller pair or exclusion.

function gcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) [a, b] = [b, a % b];
  return a || 1n;
}

function rat(n, d = 1n) {
  n = BigInt(n);
  d = BigInt(d);
  if (d < 0n) [n, d] = [-n, -d];
  const g = gcd(n, d);
  return { n: n / g, d: d / g };
}

const zero = () => rat(0);
const addR = (a, b) => rat(a.n * b.d + b.n * a.d, a.d * b.d);
const mulR = (a, b) => rat(a.n * b.n, a.d * b.d);
const negR = a => rat(-a.n, a.d);

function trim(poly) {
  while (poly.length > 1 && poly.at(-1).n === 0n) poly.pop();
  return poly;
}

function add(a, b) {
  return trim(Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
    addR(a[i] ?? zero(), b[i] ?? zero())));
}

function scale(poly, c) {
  return trim(poly.map(value => mulR(value, c)));
}

function multiply(a, b) {
  const out = Array.from({ length: a.length + b.length - 1 }, zero);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      out[i + j] = addR(out[i + j], mulR(a[i], b[j]));
    }
  }
  return trim(out);
}

function derivative(poly) {
  if (poly.length <= 1) return [zero()];
  return trim(poly.slice(1).map((value, i) => mulR(value, rat(i + 1))));
}

function zShift(poly, amount = 1) {
  return [...Array.from({ length: amount }, zero), ...poly];
}

function equal(a, b) {
  a = trim([...a]);
  b = trim([...b]);
  return a.length === b.length && a.every((x, i) =>
    x.n * b[i].d === b[i].n * x.d);
}

function show(poly) {
  return trim([...poly]).map(q => `${q.n}/${q.d}`).join(",");
}

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

// First-face data, with B(0)=0 as required by the visible U_1=0 face condition.
const A = [rat(1)];
const B = [zero(), zero(), zero(), rat(1)]; // z^3
const C = [zero(), zero(), zero(), rat(5, 3)]; // 5 z^3 / 3
const G0 = [rat(1), zero(), zero(), zero(), rat(-1)]; // 1-z^4

// The first source face is D_0 C + 5 B E_0 = 5 z^7,
// where D_0=15A-5zA' and E_0=zG_0'-5G_0.
const D0 = add(scale(A, rat(15)), scale(zShift(derivative(A)), rat(-5)));
const E0 = add(zShift(derivative(G0)), scale(G0, rat(-5)));
const firstFace = add(multiply(D0, C), scale(multiply(B, E0), rat(5)));
const firstFaceRhs = zShift([rat(5)], 7);
assert(equal(firstFace, firstFaceRhs), "first tau-face identity failed");

// For Fhat=A+tau*v*B+tau^2*v^2*F2+..., Ghat=G0+tau*v*C+tau^2*v^2*G2+...,
// choose these low-degree second-face coefficients.
const F2 = [zero(), zero(), rat(5, 3)]; // 5 z^2 / 3
const G2 = [zero(), zero(), rat(25, 9)]; // 25 z^2 / 9
const r0 = G0;
const r1 = add(C, scale(multiply(G0, B), rat(-5, 3)));
const r2 = add(
  add(add(G2, scale(multiply(G0, F2), rat(-5, 3))),
    scale(multiply(C, B), rat(-5, 3))),
  scale(multiply(G0, multiply(B, B)), rat(20, 9)));

// The tau^1 coefficient of the transformed ratio PDE, after dividing by v^2,
// is this exact polynomial identity for A=1 and the displayed first face.
const z3 = zShift([rat(1)], 3);
const z10 = zShift([rat(1)], 10);
const secondFace = add(
  add(scale(r2, rat(30)), scale(multiply(z3, r1), rat(-16))),
  add(
    multiply(B, add(r1, scale(zShift(derivative(r1)), rat(5)))),
    scale(multiply(
      multiply(zShift([rat(1)]), add(scale(F2, rat(2)), scale(multiply(B, B), rat(-1)))),
      derivative(r0)), rat(5))));
const secondFaceRhs = scale(z10, rat(-40, 3));
assert(equal(secondFace, secondFaceRhs), "second tau-face identity failed");
assert(B[0].n === 0n, "the U_1=0 face condition B(0)=0 failed");

// Equivalent axis 2-jet. Its Jacobian is 1 through the coefficient of y,
// which still does not extend it to an actual polynomial Keller pair.
const p0 = [zero(), zero(), zero(), rat(1)]; // x^3
const p1 = [rat(1)];
const p2 = [zero(), rat(5, 3)]; // 5x/3
const q0 = [zero(), rat(-1), zero(), zero(), zero(), rat(1)]; // x^5-x
const q1 = [zero(), zero(), rat(5, 3)]; // 5x^2/3
const q2 = [zero(), zero(), zero(), rat(25, 9)]; // 25x^3/9
const jacobianY0 = add(multiply(derivative(p0), q1),
  scale(multiply(p1, derivative(q0)), rat(-1)));
const jacobianY1 = add(
  add(scale(multiply(derivative(p0), q2), rat(2)),
    multiply(derivative(p1), q1)),
  add(scale(multiply(p1, derivative(q1)), rat(-1)),
    scale(multiply(p2, derivative(q0)), rat(-2))));
assert(equal(jacobianY0, [rat(1)]), "axis Jacobian constant term failed");
assert(equal(jacobianY1, [zero()]), "axis Jacobian linear term failed");

console.log("Exact distance-36 tau-face replay: PASS");
console.log(`first face: ${show(firstFace)} = ${show(firstFaceRhs)}`);
console.log(`r1: ${show(r1)}`);
console.log(`r2: ${show(r2)}`);
console.log(`second face: ${show(secondFace)} = ${show(secondFaceRhs)}`);
console.log(`axis Jacobian coefficients: y^0=${show(jacobianY0)}, y^1=${show(jacobianY1)}`);
console.log("Scope: two compatible local faces only; no all-root lift or branch exclusion.");
