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

R_in          = 36.0;   // inner arc radius of full-thickness section
R_edge        = 38.0;   // sharp leading/trailing edge arc radius
R_out         = 40.0;   // outer arc radius of full-thickness section

root_embed    = 2.0;    // all root points embedded this far inside fuselage, mm
axial_chamfer = 2.0;    // axial chamfer length at leading/trailing edge, mm

// Resolution — set PREVIEW=false for STL export.
PREVIEW  = false;

FN_BODY  = PREVIEW ? 64  : 360;  // rotate_extrude segments (circumferential)
FN_NOSE  = PREVIEW ? 32  : 360;  // ogive axial profile samples
FN_CYL   = PREVIEW ? 20  : 360;  // cylinder axial segments
                                  // 291 = ceil(570 mm / 1.96 mm) for ~1:1 quads
FN_WING  = PREVIEW ? 64  : 360;  // stabilizer arc samples

// ── Derived constants ─────────────────────────────────────────────────────────

R          = D / 2;
ogive_rho  = 8.5 * D;
total_len  = 10.0 * D;
ogive_len  = sqrt(ogive_rho^2 - (ogive_rho - R)^2);
cyl_len    = total_len - ogive_len;
root_offset = R_edge / sqrt(2);

assert(R_in   < R_edge,         "Require R_in < R_edge");
assert(R_edge < R_out,          "Require R_edge < R_out");
assert(R_out  == R,             "R_out must equal fuselage radius R");
assert(root_embed > 0,          "root_embed must be positive");
assert(2 * axial_chamfer < L,   "Require 2 * axial_chamfer < L");

echo(str("Ogive length       = ", ogive_len,         " mm"));
echo(str("Cylinder length    = ", cyl_len,           " mm"));
echo(str("Stabilizer length  = ", L,                 " mm"));
echo(str("Root embed depth   = ", root_embed,        " mm"));
echo(str("Cyl axial step     = ", cyl_len / FN_CYL,  " mm"));
echo(str("Circ step          = ", 2 * 3.14159265 * R / FN_BODY, " mm"));

// ── Fuselage ──────────────────────────────────────────────────────────────────

function ogive_r(x) =
    sqrt(ogive_rho^2 - (ogive_len - x)^2) - (ogive_rho - R);

// Profile: (r, z), tip at z=0, base at z=total_len.
//
// Three sections:
//   1. Nose tip point         [0, 0]
//   2. Ogive curve            FN_NOSE samples, uniform in x
//   3. Cylinder axial rungs   FN_CYL intermediate z-values at r=R
//      (without these the cylinder is one 570 mm tall quad → 290:1 sliver)
//   4. Base closure           [R, total_len], [0, total_len]
//
function body_profile() = concat(
    [[0, 0]],
    [for (i = [1 : FN_NOSE])
        let(x = i * ogive_len / FN_NOSE)
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

module assembly() {
    render(convexity = 10)
    rotate([0, 90, 0])
    union() {
        fuselage();
        stabilizers();
    }
}

EXPORT = "";

if      (EXPORT == "body") assembly();
else if (EXPORT == "")     assembly();
