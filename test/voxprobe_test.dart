import 'package:test/test.dart';
import 'package:voxprobe/voxprobe.dart';

void main() {
  test('ProbeSH3 directional & hemisphere irradiance baking responds strongly along light normal', () {
    final sh = ProbeSH3()
      ..addDirectionalLight(const ProbeVec3(0.0, 1.0, 0.0), const ProbeColor3(1.0, 0.8, 0.5))
      ..addAmbientHemisphere(
        skyColor: const ProbeColor3(0.2, 0.4, 0.8),
        horizonColor: const ProbeColor3(0.3, 0.3, 0.35),
        groundColor: const ProbeColor3(0.1, 0.08, 0.05),
      );

    final up = sh.evaluateIrradiance(const ProbeVec3(0.0, 1.0, 0.0));
    final down = sh.evaluateIrradiance(const ProbeVec3(0.0, -1.0, 0.0));
    final packed = sh.toFlatArray27();

    expect(up.r, greaterThan(0.8));
    expect(up.b, greaterThan(down.b * 2.0));
    expect(up.g, greaterThan(down.g * 2.0));
    expect(packed.length, 27);
  });

  test('TetrahedralProbeGrid computes normalized 5-tet barycentric weights and interpolates smoothly', () {
    final grid = TetrahedralProbeGrid.fromVolume(
      const ProbeBakeVolume(
        minCorner: ProbeVec3(-6.0, 0.0, -6.0),
        maxCorner: ProbeVec3(6.0, 6.0, 6.0),
        nx: 3,
        ny: 3,
        nz: 3,
      ),
    );
    grid.bakeFromLights([
      const ProbeLightInput.point(
        position: ProbeVec3(5.0, 3.0, 5.0),
        color: ProbeColor3(1.0, 0.3, 0.1),
        intensity: 3.0,
        range: 10.0,
      ),
    ]);

    for (final pt in const [
      ProbeVec3(-3.2, 1.4, -2.5),
      ProbeVec3(0.0, 3.0, 0.0),
      ProbeVec3(3.8, 4.1, 4.2),
      ProbeVec3(1.1, 2.2, -4.4),
    ]) {
      final loc = grid.locateSimplex(pt);
      final wSum = loc.weights.fold(0.0, (a, b) => a + b);
      expect(wSum, closeTo(1.0, 1e-6));
      for (final w in loc.weights) {
        expect(w, greaterThanOrEqualTo(-1e-6));
      }
    }

    final nearLight = grid.sampleIrradiance(const ProbeVec3(4.5, 3.0, 4.5), ProbeVec3.unitY);
    final farSide = grid.sampleIrradiance(const ProbeVec3(-5.0, 3.0, -5.0), ProbeVec3.unitY);
    expect(nearLight.r, greaterThan(farSide.r * 2.0));
  });

  test('Occluded probe filtering prevents geometry-trapped dark probes from leaking black', () {
    final grid = TetrahedralProbeGrid.fromVolume(
      const ProbeBakeVolume(
        minCorner: ProbeVec3(0.0, 0.0, 0.0),
        maxCorner: ProbeVec3(4.0, 4.0, 4.0),
        nx: 3,
        ny: 3,
        nz: 3,
      ),
    );
    grid.bakeFromLights(
      [
        const ProbeLightInput.directional(
          directionToLight: ProbeVec3(0.0, 1.0, 0.0),
          color: ProbeColor3(1.0, 1.0, 1.0),
          intensity: 1.5,
        ),
      ],
      // Occlude the exact probe at (2, 2, 2)
      isOccluded: (pos) => (pos - const ProbeVec3(2.0, 2.0, 2.0)).length < 0.25,
    );

    expect(grid.occludedProbeCount, 1);
    final sampled = grid.sampleIrradiance(const ProbeVec3(2.1, 2.1, 2.1), ProbeVec3.unitY);
    // Remaining 3 active simplex vertices preserve full brightness (~>1.0) rather than dipping!
    expect(sampled.r, greaterThan(1.0));
  });
}
