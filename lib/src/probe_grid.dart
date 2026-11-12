import 'dart:math' as math;

import 'math.dart';
import 'sh3.dart';

/// Light input used during off-line/incremental probe SH3 baking.
final class ProbeLightInput {
  final bool isDirectional;
  final ProbeVec3 positionOrDirToLight;
  final ProbeColor3 color;
  final double intensity;
  final double range;

  const ProbeLightInput.directional({
    required ProbeVec3 directionToLight,
    required this.color,
    this.intensity = 1.0,
  })  : isDirectional = true,
        positionOrDirToLight = directionToLight,
        range = double.infinity;

  const ProbeLightInput.point({
    required ProbeVec3 position,
    required this.color,
    this.intensity = 2.0,
    this.range = 12.0,
  })  : isDirectional = false,
        positionOrDirToLight = position;
}

/// A baked 3D GI probe sample storing world position, L2 SH data & validity flag.
final class LightProbeSample {
  final ProbeVec3 position;
  ProbeSH3 sh;
  bool active;

  LightProbeSample({
    required this.position,
    ProbeSH3? sh,
    this.active = true,
  }) : sh = sh ?? ProbeSH3();
}

/// Bounds and subdivision metadata for a 3D light-probe lattice volume.
final class ProbeBakeVolume {
  final ProbeVec3 minCorner;
  final ProbeVec3 maxCorner;
  final int nx;
  final int ny;
  final int nz;

  const ProbeBakeVolume({
    this.minCorner = const ProbeVec3(-10.0, 0.0, -10.0),
    this.maxCorner = const ProbeVec3(10.0, 8.0, 10.0),
    this.nx = 4,
    this.ny = 3,
    this.nz = 4,
  });
}

/// Canonical Kuhn 5-tetrahedra 3D probe volume interpolator with occlusion filtering.
final class TetrahedralProbeGrid {
  final ProbeBakeVolume volume;
  final List<LightProbeSample> probes;

  TetrahedralProbeGrid._(this.volume, this.probes);

  factory TetrahedralProbeGrid.fromVolume(ProbeBakeVolume volume) {
    final nx = math.max(2, volume.nx);
    final ny = math.max(2, volume.ny);
    final nz = math.max(2, volume.nz);
    final probes = <LightProbeSample>[];
    for (var iz = 0; iz < nz; iz++) {
      final tz = iz / (nz - 1);
      final z = volume.minCorner.z + (volume.maxCorner.z - volume.minCorner.z) * tz;
      for (var iy = 0; iy < ny; iy++) {
        final ty = iy / (ny - 1);
        final y = volume.minCorner.y + (volume.maxCorner.y - volume.minCorner.y) * ty;
        for (var ix = 0; ix < nx; ix++) {
          final tx = ix / (nx - 1);
          final x = volume.minCorner.x + (volume.maxCorner.x - volume.minCorner.x) * tx;
          probes.add(LightProbeSample(position: ProbeVec3(x, y, z)));
        }
      }
    }
    return TetrahedralProbeGrid._(
      ProbeBakeVolume(
        minCorner: volume.minCorner,
        maxCorner: volume.maxCorner,
        nx: nx,
        ny: ny,
        nz: nz,
      ),
      probes,
    );
  }

  int _idx(int ix, int iy, int iz) => (iz * volume.ny + iy) * volume.nx + ix;

  int get activeProbeCount => probes.where((p) => p.active).length;
  int get occludedProbeCount => probes.length - activeProbeCount;

  void clear() {
    for (final p in probes) {
      p.sh = ProbeSH3();
      p.active = true;
    }
  }

  /// Bakes sky hemisphere + directional/point lights into every grid probe.
  void bakeFromLights(
    List<ProbeLightInput> lights, {
    ProbeColor3 skyColor = const ProbeColor3(0.28, 0.42, 0.65),
    ProbeColor3 horizonColor = const ProbeColor3(0.45, 0.48, 0.52),
    ProbeColor3 groundColor = const ProbeColor3(0.14, 0.12, 0.1),
    bool Function(ProbeVec3 probePos)? isOccluded,
    bool Function(ProbeVec3 from, ProbeVec3 to)? isLineBlocked,
  }) {
    for (final p in probes) {
      if (isOccluded != null && isOccluded(p.position)) {
        p.active = false;
        p.sh = ProbeSH3();
        continue;
      }
      p.active = true;
      final sh = ProbeSH3()
        ..addAmbientHemisphere(
          skyColor: skyColor,
          horizonColor: horizonColor,
          groundColor: groundColor,
        );

      for (final light in lights) {
        if (light.isDirectional) {
          sh.addDirectionalLight(
            light.positionOrDirToLight,
            light.color * light.intensity,
          );
        } else {
          final toLight = light.positionOrDirToLight - p.position;
          final dist = toLight.length;
          if (dist > light.range) continue;
          if (isLineBlocked != null && isLineBlocked(p.position, light.positionOrDirToLight)) {
            continue;
          }
          final atten = math.pow(math.max(0.0, 1.0 - dist / light.range), 2.0).toDouble();
          final dir = dist > 1e-6 ? (toLight / dist) : ProbeVec3.unitY;
          sh.addDirectionalLight(dir, light.color * (light.intensity * atten));
        }
      }
      p.sh = sh;
    }
  }

  /// Computes tetrahedral simplex vertices and normalized barycentric weights at [worldPos].
  ({List<LightProbeSample> probes, List<double> weights}) locateSimplex(ProbeVec3 worldPos) {
    final nx = volume.nx;
    final ny = volume.ny;
    final nz = volume.nz;
    final sizeX = math.max(1e-6, volume.maxCorner.x - volume.minCorner.x);
    final sizeY = math.max(1e-6, volume.maxCorner.y - volume.minCorner.y);
    final sizeZ = math.max(1e-6, volume.maxCorner.z - volume.minCorner.z);

    final fx = ((worldPos.x - volume.minCorner.x) / sizeX).clamp(0.0, 1.0) * (nx - 1);
    final fy = ((worldPos.y - volume.minCorner.y) / sizeY).clamp(0.0, 1.0) * (ny - 1);
    final fz = ((worldPos.z - volume.minCorner.z) / sizeZ).clamp(0.0, 1.0) * (nz - 1);

    final ix = fx.floor().clamp(0, nx - 2);
    final iy = fy.floor().clamp(0, ny - 2);
    final iz = fz.floor().clamp(0, nz - 2);

    var u = (fx - ix).clamp(0.0, 1.0);
    var v = (fy - iy).clamp(0.0, 1.0);
    var w = (fz - iz).clamp(0.0, 1.0);

    // Alternating parity keeps adjacent cube diagonals conforming (Kuhn 5-tet)
    final flipX = (ix % 2) == 1;
    final flipY = (iy % 2) == 1;
    final flipZ = (iz % 2) == 1;
    if (flipX) u = 1.0 - u;
    if (flipY) v = 1.0 - v;
    if (flipZ) w = 1.0 - w;

    LightProbeSample corner(int cx, int cy, int cz) {
      final ox = flipX ? (1 - cx) : cx;
      final oy = flipY ? (1 - cy) : cy;
      final oz = flipZ ? (1 - cz) : cz;
      return probes[_idx(ix + ox, iy + oy, iz + oz)];
    }

    // Standard 5-tetrahedra partition of unit cube [0,1]^3:
    // Corner cuts use half-spaces u+v+w <= 1, (1-u)+(1-v)+w <= 1, etc.
    late final List<LightProbeSample> tetProbes;
    late final List<double> tetWeights;

    if (u + v + w <= 1.0) {
      // Corner at (0,0,0) with (1,0,0), (0,1,0), (0,0,1)
      tetProbes = [corner(0, 0, 0), corner(1, 0, 0), corner(0, 1, 0), corner(0, 0, 1)];
      tetWeights = [1.0 - u - v - w, u, v, w];
    } else if (u + (1.0 - v) + (1.0 - w) <= 1.0) {
      // Corner at (0,1,1) with (1,1,1), (0,0,1), (0,1,0)
      final b1 = u;
      final b2 = 1.0 - v;
      final b3 = 1.0 - w;
      tetProbes = [corner(0, 1, 1), corner(1, 1, 1), corner(0, 0, 1), corner(0, 1, 0)];
      tetWeights = [1.0 - b1 - b2 - b3, b1, b2, b3];
    } else if ((1.0 - u) + v + (1.0 - w) <= 1.0) {
      // Corner at (1,0,1) with (0,0,1), (1,1,1), (1,0,0)
      final b1 = 1.0 - u;
      final b2 = v;
      final b3 = 1.0 - w;
      tetProbes = [corner(1, 0, 1), corner(0, 0, 1), corner(1, 1, 1), corner(1, 0, 0)];
      tetWeights = [1.0 - b1 - b2 - b3, b1, b2, b3];
    } else if ((1.0 - u) + (1.0 - v) + w <= 1.0) {
      // Corner at (1,1,0) with (0,1,0), (1,0,0), (1,1,1)
      final b1 = 1.0 - u;
      final b2 = 1.0 - v;
      final b3 = w;
      tetProbes = [corner(1, 1, 0), corner(0, 1, 0), corner(1, 0, 0), corner(1, 1, 1)];
      tetWeights = [1.0 - b1 - b2 - b3, b1, b2, b3];
    } else {
      // Regular inner tetrahedron {(1,0,0), (0,1,0), (0,0,1), (1,1,1)}
      tetProbes = [corner(1, 0, 0), corner(0, 1, 0), corner(0, 0, 1), corner(1, 1, 1)];
      tetWeights = [
        0.5 * (1.0 + u - v - w),
        0.5 * (1.0 - u + v - w),
        0.5 * (1.0 - u - v + w),
        0.5 * (u + v + w - 1.0),
      ];
    }

    // Mask out inactive/occluded probes & renormalize to prevent dark light-leak artifacts
    var activeSum = 0.0;
    final filtered = List<double>.filled(4, 0.0);
    for (var i = 0; i < 4; i++) {
      if (tetProbes[i].active) {
        filtered[i] = math.max(0.0, tetWeights[i]);
        activeSum += filtered[i];
      }
    }
    if (activeSum > 1e-6) {
      for (var i = 0; i < 4; i++) {
        filtered[i] /= activeSum;
      }
    } else {
      // Fallback: nearest active probe in grid
      var bestIdx = -1;
      var bestDist2 = double.infinity;
      for (var i = 0; i < probes.length; i++) {
        if (!probes[i].active) continue;
        final d2 = probes[i].position.distanceSquaredTo(worldPos);
        if (d2 < bestDist2) {
          bestDist2 = d2;
          bestIdx = i;
        }
      }
      if (bestIdx >= 0) {
        return (probes: [probes[bestIdx]], weights: [1.0]);
      }
    }

    return (probes: tetProbes, weights: filtered);
  }

  ProbeSH3 sampleAt(ProbeVec3 worldPos) {
    final loc = locateSimplex(worldPos);
    return ProbeSH3.blend([for (final p in loc.probes) p.sh], loc.weights);
  }

  ProbeColor3 sampleIrradiance(ProbeVec3 worldPos, ProbeVec3 normal) =>
      sampleAt(worldPos).evaluateIrradiance(normal);
}
