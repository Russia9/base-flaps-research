/*
 * Parametric fuselage with deployable side arc stabilizers.
 *
 * Coordinate system after assembly:
 *   X – axial, nose at origin, positive toward the base / wake
 *   Y – lateral, first stabilizer centered on +Y
 *   Z – lateral
 *
 * All dimensions in millimetres.
 *
 */

// ── User parameters ───────────────────────────────────────────────────────────

D             = 80.0;   // fuselage outer diameter, mm
N             = 4;      // number of stabilizers
xi            = 90;     // stabilizer arc span, degrees
L             = 140.0;  // axial stabilizer chord length, mm

R_nose        = 0.5;  // spherical nose-tip radius, mm — small blunting, tangent to ogive

R_in          = 36.0;   // inner arc radius of full-thickness section
R_edge        = 38.0;   // sharp leading/trailing edge arc radius
R_out         = 40.0;   // outer arc radius of full-thickness section

root_embed    = 2.0;    // all root points embedded this far inside fuselage, mm
axial_chamfer = 2.0;    // axial chamfer length at leading/trailing edge, mm

// ── Derived constants ─────────────────────────────────────────────────────────

R          = D / 2;
ogive_rho  = 8.5 * D;
total_len  = 10.0 * D;
ogive_len  = sqrt(ogive_rho^2 - (ogive_rho - R)^2);
cyl_len    = total_len - ogive_len;
root_offset = R_edge / sqrt(2);

// Spherically blunted tangent ogive: a tip sphere of radius R_nose, internally
// tangent to the ogive arc. The ogive circle (centre (ogive_len, -(ogive_rho-R)),
// radius ogive_rho) is left untouched so it still meets the cylinder at ogive_len;
// only the sharp front is replaced by the cap, so the apex recedes to nose_apex_x.
nose_s      = sqrt((ogive_rho - R_nose)^2 - (ogive_rho - R)^2);
nose_xo     = ogive_len - nose_s;                                    // cap centre (axial)
nose_apex_x = nose_xo - R_nose;                                      // foremost point on axis
nose_tan_x  = ogive_len - ogive_rho * nose_s / (ogive_rho - R_nose); // ogive tangency (axial)
nose_tan_r  = (ogive_rho - R) * R_nose / (ogive_rho - R_nose);       // ogive tangency (radius)
nose_phi_t  = atan2(nose_tan_r, nose_tan_x - nose_xo);               // cap arc end angle, deg

assert(R_nose > 0 && R_nose < R, "Require 0 < R_nose < R");
assert(R_in   < R_edge,         "Require R_in < R_edge");
assert(R_edge < R_out,          "Require R_edge < R_out");
assert(R_out  == R,             "R_out must equal fuselage radius R");
assert(root_embed > 0,          "root_embed must be positive");
assert(2 * axial_chamfer < L,   "Require 2 * axial_chamfer < L");

echo(str("Nose tip radius    = ", R_nose,            " mm"));
echo(str("Nose apex at x      = ", nose_apex_x,       " mm"));
echo(str("Ogive length       = ", ogive_len,         " mm"));
echo(str("Cylinder length    = ", cyl_len,           " mm"));
echo(str("Stabilizer length  = ", L,                 " mm"));
echo(str("Root embed depth   = ", root_embed,        " mm"));
echo(str("Cyl axial step     = ", cyl_len / FN_CYL,  " mm"));
echo(str("Circ step          = ", 2 * 3.14159265 * R / FN_BODY, " mm"));


// Resolution — set PREVIEW=false for STL export.
PREVIEW  = false;

facet_target = 0.23;  // mm, must be < snappy L5 surface cell (~0.417mm)

FN_BODY  = PREVIEW ? 64 : ceil(2*PI*R / facet_target);
FN_CAP   = PREVIEW ? 12 : max(8, ceil(R_nose * (180 - nose_phi_t) * PI/180 / facet_target));
FN_NOSE  = PREVIEW ? 32 : ceil(ogive_len / facet_target);
FN_CYL   = PREVIEW ? 20 : ceil(cyl_len / facet_target);
FN_WING  = PREVIEW ? 64 : ceil((R_out*xi*PI/180) / facet_target);

// ── Fuselage ──────────────────────────────────────────────────────────────────

function ogive_r(x) =
    sqrt(ogive_rho^2 - (ogive_len - x)^2) - (ogive_rho - R);

// Profile: (r, z), blunt apex at z=nose_apex_x, base at z=total_len.
//
// Four sections:
//   1. Spherical nose cap     FN_CAP samples, apex on axis → ogive tangency
//   2. Ogive curve            FN_NOSE samples, tangency → base, uniform in x
//   3. Cylinder axial rungs   FN_CYL intermediate z-values at r=R
//      (without these the cylinder is one 570 mm tall quad → 290:1 sliver)
//   4. Base closure           [R, total_len], [0, total_len]
//
function body_profile() = concat(
    [for (i = [0 : FN_CAP])
        let(phi = 180 - (180 - nose_phi_t) * i / FN_CAP)
        [R_nose * sin(phi), nose_xo + R_nose * cos(phi)]
    ],
    [for (i = [1 : FN_NOSE])
        let(x = nose_tan_x + i * (ogive_len - nose_tan_x) / FN_NOSE)
        [ogive_r(x), x]
    ],
    [for (j = [1 : FN_CYL - 1])
        [R, ogive_len + j * cyl_len / FN_CYL]
    ],
    [[R, total_len],
     [0, total_len]]
);

module fuselage() {
    rotate_extrude($fn = FN_BODY)
        polygon(body_profile());
}

// ── Stabilizer geometry ───────────────────────────────────────────────────────

wing_center  = [-root_offset, R_out + root_offset];
root_y       = R - root_embed;   // all root points embedded root_embed below surface

function pt_on_arc(r, a) =
    [wing_center[0] + r * cos(a), wing_center[1] + r * sin(a)];

function angle_of(p) =
    atan2(p[1] - wing_center[1], p[0] - wing_center[0]);

function x_at_y(r, y) =
    wing_center[0] + sqrt(max(0, r^2 - (y - wing_center[1])^2));

root_inner = [x_at_y(R_in,   root_y), root_y];
root_outer = [x_at_y(R_out,  root_y), root_y];
root_edge  = [x_at_y(R_edge, root_y), root_y];   // embedded, not flush

theta_inner = angle_of(root_inner);
theta_outer = angle_of(root_outer);
theta_edge  = angle_of(root_edge);
theta_tip   = theta_edge + xi;

function edge_pt(i)  = pt_on_arc(R_edge, theta_edge  + (theta_tip - theta_edge)  * i / FN_WING);
function outer_pt(i) = pt_on_arc(R_out,  theta_outer + (theta_tip - theta_outer) * i / FN_WING);
function inner_pt(i) = pt_on_arc(R_in,   theta_inner + (theta_tip - theta_inner) * i / FN_WING);

function p3(p, z) = [p[0], p[1], z];
function vi(block, i) = block * (FN_WING + 1) + i;

// ── Stabilizer solid ──────────────────────────────────────────────────────────
//
// Blocks 0–5, each (FN_WING+1) points:
//   0  z0  edge arc  (leading knife-edge)
//   1  z1  outer arc
//   2  z1  inner arc
//   3  z2  outer arc
//   4  z2  inner arc
//   5  z3  edge arc  (trailing knife-edge)
//
// Winding: CCW from outside (right-hand outward normal).

module stabilizer_one() {
    z0 = total_len - L;
    z1 = z0 + axial_chamfer;
    z2 = total_len - axial_chamfer;
    z3 = total_len;
    n  = FN_WING;

    pts = concat(
        [for (i=[0:n]) p3(edge_pt(i),  z0)],
        [for (i=[0:n]) p3(outer_pt(i), z1)],
        [for (i=[0:n]) p3(inner_pt(i), z1)],
        [for (i=[0:n]) p3(outer_pt(i), z2)],
        [for (i=[0:n]) p3(inner_pt(i), z2)],
        [for (i=[0:n]) p3(edge_pt(i),  z3)]
    );

    outer_lead   = [for (i=[0:n-1]) each [[vi(0,i),vi(1,i),vi(1,i+1)],   [vi(0,i),vi(1,i+1),vi(0,i+1)]]];
    inner_lead   = [for (i=[0:n-1]) each [[vi(0,i+1),vi(2,i+1),vi(2,i)], [vi(0,i+1),vi(2,i),vi(0,i)]]];
    outer_barrel = [for (i=[0:n-1]) each [[vi(1,i),vi(3,i),vi(3,i+1)],   [vi(1,i),vi(3,i+1),vi(1,i+1)]]];
    inner_barrel = [for (i=[0:n-1]) each [[vi(2,i+1),vi(4,i+1),vi(4,i)], [vi(2,i+1),vi(4,i),vi(2,i)]]];
    outer_trail  = [for (i=[0:n-1]) each [[vi(3,i),vi(5,i),vi(5,i+1)],   [vi(3,i),vi(5,i+1),vi(3,i+1)]]];
    inner_trail  = [for (i=[0:n-1]) each [[vi(4,i+1),vi(5,i+1),vi(5,i)], [vi(4,i+1),vi(5,i),vi(4,i)]]];

    root_cap = [
        [vi(0,0), vi(2,0), vi(1,0)],
        [vi(1,0), vi(2,0), vi(4,0)],
        [vi(1,0), vi(4,0), vi(3,0)],
        [vi(5,0), vi(3,0), vi(4,0)]
    ];

    tip_cap = [
        [vi(0,n), vi(1,n), vi(2,n)],
        [vi(1,n), vi(3,n), vi(4,n)],
        [vi(1,n), vi(4,n), vi(2,n)],
        [vi(5,n), vi(4,n), vi(3,n)]
    ];

    polyhedron(
        points   = pts,
        faces    = concat(
            outer_lead, inner_lead,
            outer_barrel, inner_barrel,
            outer_trail, inner_trail,
            root_cap, tip_cap
        ),
        convexity = 10
    );
}

module stabilizers() {
    for (k = [0 : N - 1])
        rotate([0, 0, k * 360 / N])
            stabilizer_one();
}

// ── Assembly ──────────────────────────────────────────────────────────────────

EXPORT = "";

module assembly(part = "all") {
    render(convexity = 10)
    rotate([0, 90, 0])
    union() {
        if (part == "all" || part == "fuselage")    fuselage();
        if (part == "all" || part == "stabilizers") stabilizers();
    }
}

if      (EXPORT == "fuselage")    assembly("fuselage");
else if (EXPORT == "stabilizers") assembly("stabilizers");
else                              assembly("all");
