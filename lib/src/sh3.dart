import 'dart:math' as math;
import 'dart:typed_data';

import 'math.dart';

/// 9-band L2 RGB Spherical Harmonics container for real-time diffuse irradiance evaluation.
final class ProbeSH3 {
  final List<ProbeColor3> coeffs;

  ProbeSH3([List<ProbeColor3>? initial])
      : coeffs = initial != null && initial.length == 9
            ? List<ProbeColor3>.of(initial, growable: false)
            : List<ProbeColor3>.filled(9, ProbeColor3.black, growable: false);

  static List<double> _basis(ProbeVec3 dir) {
    final d = dir.normalized;
    final x = d.x;
    final y = d.y;
    final z = d.z;
    return [
      0.282095, // L00
      0.488603 * y, // L1,-1 (Y)
      0.488603 * z, // L1,0  (Z)
      0.488603 * x, // L1,1  (X)
      1.092548 * x * y, // L2,-2
      1.092548 * y * z, // L2,-1
      0.315392 * (3.0 * z * z - 1.0), // L2,0
      1.092548 * x * z, // L2,1
      0.546274 * (x * x - y * y), // L2,2
    ];
  }

  /// Projects a directional irradiance impulse arriving from [directionToLight] into L2 SH.
  void addDirectionalLight(ProbeVec3 directionToLight, ProbeColor3 irradiance) {
    final y = _basis(directionToLight);
    // Ramamoorthi-Hanrahan cosine convolution multipliers A0 = pi, A1 = 2pi/3, A2 = pi/4
    const a0 = math.pi;
    const a1 = 2.0 * math.pi / 3.0;
    const a2 = math.pi * 0.25;
    const bandA = [a0, a1, a1, a1, a2, a2, a2, a2, a2];
    for (var i = 0; i < 9; i++) {
      coeffs[i] = coeffs[i] + irradiance * (y[i] * bandA[i]);
    }
  }

  /// Adds an analytical three-band environment sky/horizon/ground gradient.
  void addAmbientHemisphere({
    required ProbeColor3 skyColor,
    required ProbeColor3 horizonColor,
    required ProbeColor3 groundColor,
  }) {
    final l00 = (skyColor + (horizonColor * 2.0) + groundColor) * (0.282095 * 0.25 * math.pi);
    final l1y = (skyColor - groundColor) * (0.488603 * 0.333333 * math.pi);
    final l2y = ((skyColor + groundColor) - (horizonColor * 2.0)) * (0.315392 * 0.25 * math.pi);
    // Project using matching basis slots: L00 + L1Y (index 1) + L2 zonal aligned along +Y
    coeffs[0] = coeffs[0] + l00;
    coeffs[1] = coeffs[1] + l1y;
    // Decompose (3y^2 - 1) into {L20(3z^2-1), L22(x^2-y^2)}:
    // 3y^2 - 1 = -0.5*(3z^2 - 1) - 1.5*(x^2 - y^2) on x^2+y^2+z^2=1
    coeffs[6] = coeffs[6] + l2y * (-0.5);
    coeffs[8] = coeffs[8] + l2y * (-0.5 * 1.7320508);
  }

  /// Evaluates smooth clamped diffuse irradiance for a surface normal [normal].
  ProbeColor3 evaluateIrradiance(ProbeVec3 normal) {
    final y = _basis(normal);
    var sum = ProbeColor3.black;
    for (var i = 0; i < 9; i++) {
      sum = sum + coeffs[i] * y[i];
    }
    return sum.clamped(0.0, 1e6);
  }

  ProbeSH3 operator +(ProbeSH3 other) => ProbeSH3([
        for (var i = 0; i < 9; i++) coeffs[i] + other.coeffs[i],
      ]);

  ProbeSH3 operator *(double scalar) => ProbeSH3([
        for (var i = 0; i < 9; i++) coeffs[i] * scalar,
      ]);

  static ProbeSH3 blend(List<ProbeSH3> samples, List<double> weights) {
    final out = ProbeSH3();
    for (var k = 0; k < samples.length && k < weights.length; k++) {
      final w = weights[k];
      if (w.abs() < 1e-9) continue;
      final sc = samples[k].coeffs;
      for (var i = 0; i < 9; i++) {
        out.coeffs[i] = out.coeffs[i] + sc[i] * w;
      }
    }
    return out;
  }

  /// Packs all 9 RGB coefficients into a 27-float buffer for shader uniforms.
  Float32List toFlatArray27() {
    final out = Float32List(27);
    for (var i = 0; i < 9; i++) {
      final c = coeffs[i];
      out[i * 3] = c.r;
      out[i * 3 + 1] = c.g;
      out[i * 3 + 2] = c.b;
    }
    return out;
  }
}
