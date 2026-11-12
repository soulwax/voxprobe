import 'dart:math' as math;

/// Lightweight 3D position/direction vector for probe math.
final class ProbeVec3 {
  final double x;
  final double y;
  final double z;

  const ProbeVec3(this.x, this.y, this.z);

  static const zero = ProbeVec3(0.0, 0.0, 0.0);
  static const unitY = ProbeVec3(0.0, 1.0, 0.0);

  ProbeVec3 operator +(ProbeVec3 o) => ProbeVec3(x + o.x, y + o.y, z + o.z);
  ProbeVec3 operator -(ProbeVec3 o) => ProbeVec3(x - o.x, y - o.y, z - o.z);
  ProbeVec3 operator *(double s) => ProbeVec3(x * s, y * s, z * s);
  ProbeVec3 operator /(double s) => ProbeVec3(x / s, y / s, z / s);

  double dot(ProbeVec3 o) => x * o.x + y * o.y + z * o.z;

  ProbeVec3 cross(ProbeVec3 o) => ProbeVec3(
        y * o.z - z * o.y,
        z * o.x - x * o.z,
        x * o.y - y * o.x,
      );

  double get lengthSquared => x * x + y * y + z * z;
  double get length => math.sqrt(lengthSquared);

  ProbeVec3 get normalized {
    final len = length;
    if (len < 1e-9) return unitY;
    return this / len;
  }

  double distanceSquaredTo(ProbeVec3 o) => (this - o).lengthSquared;
  double distanceTo(ProbeVec3 o) => (this - o).length;
}

/// Linear RGB irradiance triple.
final class ProbeColor3 {
  final double r;
  final double g;
  final double b;

  const ProbeColor3(this.r, this.g, this.b);

  static const black = ProbeColor3(0.0, 0.0, 0.0);
  static const white = ProbeColor3(1.0, 1.0, 1.0);

  ProbeColor3 operator +(ProbeColor3 o) => ProbeColor3(r + o.r, g + o.g, b + o.b);
  ProbeColor3 operator -(ProbeColor3 o) => ProbeColor3(r - o.r, g - o.g, b - o.b);
  ProbeColor3 operator *(double s) => ProbeColor3(r * s, g * s, b * s);

  ProbeColor3 clamped([double min = 0.0, double max = 1e6]) => ProbeColor3(
        r.clamp(min, max),
        g.clamp(min, max),
        b.clamp(min, max),
      );

  static ProbeColor3 lerp(ProbeColor3 a, ProbeColor3 b, double t) =>
      ProbeColor3(a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t);
}
