/*
 * Parametric fuselage with deployable side arc stabilizers.
 *
 * Coordinate system after assembly:
 *   X - axial, nose at origin, positive toward the base and wake
 *   Y - lateral, first stabilizer centered on +Y
 *   Z - lateral
 *
 * Dimensions are in millimeters. The fuselage is the same ogive-cylinder used
 * by geometry/model.scad. Stabilizers are cylindrical shell panels hinged to
 * the fuselage and swung outward into the flow.
 */

// Parameters overridable from the OpenSCAD CLI.
D            = 80.0;   // fuselage diameter
N            = 4;      // number of stabilizers
xi           = 45;     // stabilizer arc angle, degrees
L            = 140.0;  // axial stabilizer length
R_in         = 36.0;   // inner radius of full-thickness section
R_edge       = 38.0;   // radius of sharp leading/trailing edge
R_out        = 40.0;   // outer radius of full-thickness section
chamfer_len  = 2.0;    // 2 mm axial taper gives 45 deg chamfers
root_overlap = 1.0;    // finite overlap into fuselage avoids line-contact STL
hinge_radius = R_out - root_overlap;
deploy_angle = 90.0;   // outward swing angle from the stowed shell position

// Derived fuselage constants.
R          = D / 2;
ogive_rho  = 8.5 * D;
total_len  = 10.0 * D;
ogive_len  = sqrt(ogive_rho * ogive_rho - (ogive_rho - R) * (ogive_rho - R));
cyl_len    = total_len - ogive_len;

echo(str("Ogive length       = ", ogive_len, " mm"));
echo(str("Cylinder length    = ", cyl_len, " mm"));
echo(str("Stabilizer length  = ", L, " mm"));
echo(str("Deploy angle       = ", deploy_angle, " degrees"));

// Resolution controls.
FN_BODY = 360;
FN_NOSE = 360;
FN_WING = 360;

function ogive_r(x) =
    sqrt(ogive_rho * ogive_rho - (ogive_len - x) * (ogive_len - x))
    - (ogive_rho - R);

function body_profile() = concat(
    [[0, 0]],
    [for (i = [1 : FN_NOSE])
        let(x = i * ogive_len / FN_NOSE)
        [ogive_r(x), x]
    ],
    [[R, total_len],
     [0, total_len]]
);

// Z-axis-aligned body; tip at Z=0, base at Z=total_len.
module fuselage() {
    rotate_extrude($fn = FN_BODY)
    polygon(body_profile());
}

wing_start = total_len - L;
wing_end   = total_len;

// Radius-vs-axial profile for the shell panel before deployment. Each axial
// end collapses to R_edge, forming a sharp leading/trailing edge instead of a
// blunt surface. With R36/R38/R40 and chamfer_len=2 mm, the chamfer faces are
// 45 degrees. The whole panel is swung outward about a long hinge axis slightly
// inside the fuselage, giving a finite solid overlap at the root. The visible
// root remains at the OY/OZ fuselage intersections after the boolean union.
module stabilizer_profile() {
    assert(R_in < R_edge && R_edge < R_out,
        "Require R_in < R_edge < R_out");
    assert(2 * chamfer_len < L,
        "Require 2 * chamfer_len < L");

    polygon([
        [R_edge, wing_start],
        [R_out,  wing_start + chamfer_len],
        [R_out,  wing_end - chamfer_len],
        [R_edge, wing_end],
        [R_in,   wing_end - chamfer_len],
        [R_in,   wing_start + chamfer_len]
    ]);
}

module rotate_about_z_at(angle, point) {
    translate(point)
    rotate([0, 0, angle])
    translate([-point[0], -point[1], -point[2]])
    children();
}

module stowed_stabilizer_one() {
    rotate([0, 0, 90])
    rotate_extrude(angle = xi, $fn = FN_WING)
    stabilizer_profile();
}

module stabilizer_one() {
    hinge_angle = 90;
    hinge_point = [
        hinge_radius * cos(hinge_angle),
        hinge_radius * sin(hinge_angle),
        0
    ];

    rotate_about_z_at(-deploy_angle, hinge_point)
    stowed_stabilizer_one();
}

module stabilizers() {
    for (k = [0 : N - 1])
        rotate([0, 0, k * 360 / N])
        stabilizer_one();
}

module assembly() {
    render(convexity = 10)
    rotate([0, 90, 0])
    union() {
        fuselage();
        stabilizers();
    }
}

EXPORT = "";

module body() { assembly(); }

if      (EXPORT == "body") body();
else if (EXPORT == "")     body();
