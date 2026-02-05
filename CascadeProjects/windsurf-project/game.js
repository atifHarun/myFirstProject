// Game variables
const canvas = document.getElementById('gameCanvas');
const ctx = canvas.getContext('2d');
const startBtn = document.getElementById('startBtn');
const pauseBtn = document.getElementById('pauseBtn');
const resetBtn = document.getElementById('resetBtn');
const scoreValue = document.getElementById('scoreValue');

// Game state
let gameRunning = false;
let gamePaused = false;
let score = 0;
let animationId = null;
let gameOver = false;

// Physics constants
const GRAVITY = 0.3;
const FRICTION = 0.99;
const BOUNCE_DAMPING = 0.7;

const SIZE_MERGE_EPSILON = 0.1;

const SHAPES = [
    { key: 'CIRCLE', name: 'Circle', sides: null, size: 16, color: '#FF6B6B' },
    { key: 'TRIANGLE', name: 'Triangle', sides: 3, size: 20, color: '#F7B801' },
    { key: 'QUADRILATERAL', name: 'Quadrilateral', sides: 4, size: 24, color: '#4ECDC4' },
    { key: 'PENTAGON', name: 'Pentagon', sides: 5, size: 28, color: '#45B7D1' },
    { key: 'HEXAGON', name: 'Hexagon', sides: 6, size: 32, color: '#6A4C93' },
    { key: 'HEPTAGON', name: 'Heptagon', sides: 7, size: 36, color: '#FF9F1C' },
    { key: 'OCTAGON', name: 'Octagon', sides: 8, size: 40, color: '#2EC4B6' },
    { key: 'NONAGON', name: 'Nonagon', sides: 9, size: 44, color: '#E71D36' },
    { key: 'DECAGON', name: 'Decagon', sides: 10, size: 48, color: '#3A86FF' },
    { key: 'HENDECAGON', name: 'Hendecagon', sides: 11, size: 52, color: '#8338EC' },
    { key: 'DODECAGON', name: 'Dodecagon', sides: 12, size: 56, color: '#06D6A0' },
    { key: 'TRISKAIDECAGON', name: 'Triskaidecagon', sides: 13, size: 60, color: '#EF476F' },
    { key: 'TETRADECAGON', name: 'Tetradecagon', sides: 14, size: 64, color: '#118AB2' }
];

function clampNearZero(obj) {
    if (Math.abs(obj.vx) < 0.02) obj.vx = 0;
    if (Math.abs(obj.vy) < 0.02) obj.vy = 0;
    if (Math.abs(obj.angularVelocity) < 0.01) obj.angularVelocity = 0;
}

const SHAPE_ORDER = SHAPES.map((s) => s.key);
const SHAPE_CONFIG = Object.fromEntries(SHAPES.map((s) => [s.key, s]));
const POLY_SIDES = Object.fromEntries(SHAPES.filter((s) => s.sides).map((s) => [s.key, s.sides]));
const SHAPE_LEVEL = Object.fromEntries(SHAPE_ORDER.map((key, idx) => [key, idx + 1]));
const OBJECT_RESTITUTION = 0.08;
const WALL_RESTITUTION = 0.05;
const FLOOR_RESTITUTION = 0.03;
const STOP_BOUNCE_THRESHOLD = 0.6;
const ANGULAR_DAMPING = 0.0;
const SURFACE_FRICTION = 0.15;
const ANGULAR_SPEED_LIMIT = 0.0;

// Object types with different colors and sizes
const OBJECT_TYPES = Object.fromEntries(
    SHAPES.map((s) => [s.key, { color: s.color, name: s.name }])
);

// Game objects array
let gameObjects = [];

// Player-controlled held object (preview at the top)
let heldObject = null;
let heldX = 0;

const SPAWN_TYPES = ['CIRCLE', 'TRIANGLE', 'QUADRILATERAL', 'PENTAGON'];
const TOP_GAME_OVER_BOUNDARY = 80;
const REST_VELOCITY_THRESHOLD = 0.22;
const REST_ANGULAR_THRESHOLD = 0.05;
const REST_FRAMES_REQUIRED = 12;
const SLEEP_FRAMES_REQUIRED = 22;
const WAKE_VELOCITY_THRESHOLD = 0.5;
const SOLVER_ITERATIONS = 6;
const PENETRATION_SLOP = 0.5;
const POSITION_CORRECTION_PERCENT = 0.7;
const FLOOR_SETTLE_DAMPING = 0.7;

function drawTopBoundaryLine() {
    ctx.save();
    ctx.strokeStyle = 'rgba(0, 0, 0, 0.6)';
    ctx.lineWidth = 2;
    ctx.setLineDash([6, 6]);
    ctx.beginPath();
    ctx.moveTo(0, TOP_GAME_OVER_BOUNDARY);
    ctx.lineTo(canvas.width, TOP_GAME_OVER_BOUNDARY);
    ctx.stroke();
    ctx.restore();
}

// Tracks which circle pairs were colliding last frame (to avoid spamming the console)
let activeCircleCollisions = new Set();

function getCirclePairKey(a, b) {
    // Stable key regardless of order
    const ax = a.x.toFixed(1);
    const ay = a.y.toFixed(1);
    const bx = b.x.toFixed(1);
    const by = b.y.toFixed(1);
    const key1 = `${ax},${ay}|${bx},${by}`;
    const key2 = `${bx},${by}|${ax},${ay}`;
    return key1 < key2 ? key1 : key2;
}

function handleCircleCollisions() {
    const newActive = new Set();

    // Reset collision state and color each frame
    for (const obj of gameObjects) {
        obj.isColliding = false;
        obj.color = obj.baseColor;
    }

    for (let i = 0; i < gameObjects.length; i++) {
        const a = gameObjects[i];
        if (a.type !== 'CIRCLE') continue;

        for (let j = i + 1; j < gameObjects.length; j++) {
            const b = gameObjects[j];
            if (b.type !== 'CIRCLE') continue;

            if (a.collidesWith(b)) {
                a.isColliding = true;
                b.isColliding = true;
                a.color = '#FFD93D';
                b.color = '#FFD93D';

                const key = getCirclePairKey(a, b);
                newActive.add(key);
                if (!activeCircleCollisions.has(key)) {
                    console.log('Circle collision detected!', {
                        a: { x: a.x, y: a.y, r: a.size },
                        b: { x: b.x, y: b.y, r: b.size }
                    });
                }
            }
        }
    }

    activeCircleCollisions = newActive;
}

function handleSameSizeMerges() {
    // Find all pairs that should merge, then apply merges after scanning.
    // This avoids modifying the array while iterating over it.
    const indicesToRemove = new Set();
    const objectsToAdd = [];

    for (let i = 0; i < gameObjects.length; i++) {
        if (indicesToRemove.has(i)) continue;
        const a = gameObjects[i];

        for (let j = i + 1; j < gameObjects.length; j++) {
            if (indicesToRemove.has(j)) continue;
            const b = gameObjects[j];

            if (a.type !== b.type) continue;

            const sameSize = Math.abs(a.size - b.size) <= SIZE_MERGE_EPSILON;
            if (!sameSize) continue;

            if (!detectCollision(a, b)) continue;

            const currentIndex = SHAPE_ORDER.indexOf(a.type);
            if (currentIndex === -1) continue;
            if (currentIndex >= SHAPE_ORDER.length - 1) continue;

            const nextType = SHAPE_ORDER[currentIndex + 1];
            const newSize = SHAPE_CONFIG[nextType]?.size ?? a.size;

            const level = SHAPE_LEVEL[a.type] ?? 0;
            score += level + level;

            indicesToRemove.add(i);
            indicesToRemove.add(j);

            const x = (a.x + b.x) / 2;
            const y = (a.y + b.y) / 2;

            const merged = new GameObject(x, y, nextType, newSize);
            merged.vx = (a.vx + b.vx) / 2;
            // Give merged objects a little "pop" so forced merges feel bouncy
            merged.vy = Math.min((a.vy + b.vy) / 2, 2) - 3;
            merged.angularVelocity = 0;
            objectsToAdd.push(merged);

            // One merge per object per frame
            break;
        }
    }

    if (indicesToRemove.size === 0) return;

    gameObjects = gameObjects.filter((_, idx) => !indicesToRemove.has(idx));
    gameObjects.push(...objectsToAdd);
}

function rotatePoint(x, y, angle) {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return { x: x * c - y * s, y: x * s + y * c };
}

function getPolygonVertices(obj) {
    const sides = POLY_SIDES[obj.type];
    if (!sides) return null;

    const r = obj.size;
    const local = [];
    for (let i = 0; i < sides; i++) {
        const a = (Math.PI * 2 * i) / sides - Math.PI / 2;
        local.push({ x: Math.cos(a) * r, y: Math.sin(a) * r });
    }

    return local.map((p) => {
        const rp = rotatePoint(p.x, p.y, obj.angle);
        return { x: rp.x + obj.x, y: rp.y + obj.y };
    });
}

function normalize(x, y) {
    const len = Math.sqrt(x * x + y * y);
    if (len === 0) return { x: 1, y: 0 };
    return { x: x / len, y: y / len };
}

function dot(ax, ay, bx, by) {
    return ax * bx + ay * by;
}

function projectPolygon(vertices, axis) {
    let min = dot(vertices[0].x, vertices[0].y, axis.x, axis.y);
    let max = min;
    for (let i = 1; i < vertices.length; i++) {
        const p = dot(vertices[i].x, vertices[i].y, axis.x, axis.y);
        if (p < min) min = p;
        if (p > max) max = p;
    }
    return { min, max };
}

function projectCircle(circle, axis) {
    const center = dot(circle.x, circle.y, axis.x, axis.y);
    return { min: center - circle.size, max: center + circle.size };
}

function intervalOverlap(a, b) {
    const overlap = Math.min(a.max, b.max) - Math.max(a.min, b.min);
    return overlap;
}

function getPolygonAxes(vertices) {
    const axes = [];
    for (let i = 0; i < vertices.length; i++) {
        const p1 = vertices[i];
        const p2 = vertices[(i + 1) % vertices.length];
        const edgeX = p2.x - p1.x;
        const edgeY = p2.y - p1.y;
        axes.push(normalize(-edgeY, edgeX));
    }
    return axes;
}

function closestPointOnSegment(px, py, ax, ay, bx, by) {
    const abx = bx - ax;
    const aby = by - ay;
    const apx = px - ax;
    const apy = py - ay;
    const abLen2 = abx * abx + aby * aby;
    const t = abLen2 === 0 ? 0 : Math.max(0, Math.min(1, (apx * abx + apy * aby) / abLen2));
    return { x: ax + abx * t, y: ay + aby * t };
}

function closestPointOnPolygon(circle, vertices) {
    let best = null;
    let bestDist2 = Infinity;
    for (let i = 0; i < vertices.length; i++) {
        const a = vertices[i];
        const b = vertices[(i + 1) % vertices.length];
        const p = closestPointOnSegment(circle.x, circle.y, a.x, a.y, b.x, b.y);
        const dx = circle.x - p.x;
        const dy = circle.y - p.y;
        const d2 = dx * dx + dy * dy;
        if (d2 < bestDist2) {
            bestDist2 = d2;
            best = p;
        }
    }
    return best;
}

function supportPoint(vertices, dirX, dirY) {
    let best = vertices[0];
    let bestDot = dot(best.x, best.y, dirX, dirY);
    for (let i = 1; i < vertices.length; i++) {
        const v = vertices[i];
        const d = dot(v.x, v.y, dirX, dirY);
        if (d > bestDot) {
            bestDot = d;
            best = v;
        }
    }
    return best;
}

function detectCollision(a, b) {
    const aPoly = getPolygonVertices(a);
    const bPoly = getPolygonVertices(b);

    // Circle vs Circle
    if (!aPoly && !bPoly) {
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        const dist2 = dx * dx + dy * dy;
        const minDist = a.size + b.size;
        if (dist2 >= minDist * minDist) return null;
        const dist = Math.sqrt(dist2) || 1;
        const normal = { x: dx / dist, y: dy / dist };
        return {
            normal,
            penetration: minDist - dist,
            contact: { x: a.x + normal.x * a.size, y: a.y + normal.y * a.size }
        };
    }

    // Polygon vs Polygon
    if (aPoly && bPoly) {
        const axes = [...getPolygonAxes(aPoly), ...getPolygonAxes(bPoly)];
        let smallestOverlap = Infinity;
        let smallestAxis = null;

        for (const axis of axes) {
            const pA = projectPolygon(aPoly, axis);
            const pB = projectPolygon(bPoly, axis);
            const o = intervalOverlap(pA, pB);
            if (o <= 0) return null;
            if (o < smallestOverlap) {
                smallestOverlap = o;
                smallestAxis = axis;
            }
        }

        const centerDir = normalize(b.x - a.x, b.y - a.y);
        if (dot(smallestAxis.x, smallestAxis.y, centerDir.x, centerDir.y) < 0) {
            smallestAxis = { x: -smallestAxis.x, y: -smallestAxis.y };
        }

        const contact = supportPoint(aPoly, smallestAxis.x, smallestAxis.y);
        return { normal: smallestAxis, penetration: smallestOverlap, contact };
    }

    // Circle vs Polygon
    const circle = aPoly ? b : a;
    const polyObj = aPoly ? a : b;
    const poly = aPoly ? aPoly : bPoly;
    const axes = getPolygonAxes(poly);
    const closest = closestPointOnPolygon(circle, poly);
    if (closest) {
        axes.push(normalize(circle.x - closest.x, circle.y - closest.y));
    }

    let smallestOverlap = Infinity;
    let smallestAxis = null;

    for (const axis of axes) {
        const pPoly = projectPolygon(poly, axis);
        const pCircle = projectCircle(circle, axis);
        const o = intervalOverlap(pPoly, pCircle);
        if (o <= 0) return null;
        if (o < smallestOverlap) {
            smallestOverlap = o;
            smallestAxis = axis;
        }
    }

    const dir = normalize(polyObj.x - circle.x, polyObj.y - circle.y);
    // smallestAxis currently points from poly to circle or circle to poly depending on construction
    // Make it point from a to b
    let normal = smallestAxis;
    if (aPoly) {
        // a is polygon, b is circle -> normal should point from a to b
        const cd = normalize(circle.x - polyObj.x, circle.y - polyObj.y);
        if (dot(normal.x, normal.y, cd.x, cd.y) < 0) normal = { x: -normal.x, y: -normal.y };
    } else {
        // a is circle, b is polygon -> normal should point from a to b
        const cd = normalize(polyObj.x - circle.x, polyObj.y - circle.y);
        if (dot(normal.x, normal.y, cd.x, cd.y) < 0) normal = { x: -normal.x, y: -normal.y };
    }

    const contact = aPoly
        ? { x: circle.x - normal.x * circle.size, y: circle.y - normal.y * circle.size }
        : { x: circle.x + normal.x * circle.size, y: circle.y + normal.y * circle.size };

    return { normal, penetration: smallestOverlap, contact };
}

function resolveObjectCollisions() {
    for (let i = 0; i < gameObjects.length; i++) {
        const a = gameObjects[i];

        for (let j = i + 1; j < gameObjects.length; j++) {
            const b = gameObjects[j];
            const hit = detectCollision(a, b);
            if (!hit) continue;

            const nx = hit.normal.x;
            const ny = hit.normal.y;

            const invMassA = a.isSleeping ? 0 : 1 / a.mass;
            const invMassB = b.isSleeping ? 0 : 1 / b.mass;
            const invMassSum = invMassA + invMassB;
            if (invMassSum === 0) continue;

            // Positional correction
            const correction =
                (Math.max(0, hit.penetration - PENETRATION_SLOP) * POSITION_CORRECTION_PERCENT) / invMassSum;
            a.x -= nx * correction * invMassA;
            a.y -= ny * correction * invMassA;
            b.x += nx * correction * invMassB;
            b.y += ny * correction * invMassB;

            const cx = hit.contact.x;
            const cy = hit.contact.y;
            const raX = cx - a.x;
            const raY = cy - a.y;
            const rbX = cx - b.x;
            const rbY = cy - b.y;

            const velAx = a.vx + (-a.angularVelocity * raY);
            const velAy = a.vy + (a.angularVelocity * raX);
            const velBx = b.vx + (-b.angularVelocity * rbY);
            const velBy = b.vy + (b.angularVelocity * rbX);

            const rvx = velBx - velAx;
            const rvy = velBy - velAy;

            const velAlongNormal = rvx * nx + rvy * ny;
            if (velAlongNormal > 0) continue;

            const raCrossN = raX * ny - raY * nx;
            const rbCrossN = rbX * ny - rbY * nx;
            const invInertiaA = a.isSleeping ? 0 : 1 / a.inertia;
            const invInertiaB = b.isSleeping ? 0 : 1 / b.inertia;

            const denom = invMassSum + (raCrossN * raCrossN) * invInertiaA + (rbCrossN * rbCrossN) * invInertiaB;
            const jn = denom === 0 ? 0 : (-(1 + OBJECT_RESTITUTION) * velAlongNormal) / denom;

            const impX = jn * nx;
            const impY = jn * ny;

            if (!a.isSleeping) {
                a.vx -= impX * invMassA;
                a.vy -= impY * invMassA;
            }
            if (!b.isSleeping) {
                b.vx += impX * invMassB;
                b.vy += impY * invMassB;
            }

            if (!a.isSleeping) a.angularVelocity -= raCrossN * jn * invInertiaA;
            if (!b.isSleeping) b.angularVelocity += rbCrossN * jn * invInertiaB;

            wakeIfNeeded(a);
            wakeIfNeeded(b);

            const tx = -ny;
            const ty = nx;
            const velAlongTangent = rvx * tx + rvy * ty;

            const raCrossT = raX * ty - raY * tx;
            const rbCrossT = rbX * ty - rbY * tx;
            const denomT = invMassSum + (raCrossT * raCrossT) * invInertiaA + (rbCrossT * rbCrossT) * invInertiaB;
            let jt = denomT === 0 ? 0 : (-velAlongTangent) / denomT;

            const maxFriction = Math.abs(jn) * SURFACE_FRICTION;
            if (jt > maxFriction) jt = maxFriction;
            if (jt < -maxFriction) jt = -maxFriction;

            const fX = jt * tx;
            const fY = jt * ty;

            if (!a.isSleeping) {
                a.vx -= fX * invMassA;
                a.vy -= fY * invMassA;
            }
            if (!b.isSleeping) {
                b.vx += fX * invMassB;
                b.vy += fY * invMassB;
            }

            if (!a.isSleeping) a.angularVelocity -= raCrossT * jt * invInertiaA;
            if (!b.isSleeping) b.angularVelocity += rbCrossT * jt * invInertiaB;

            clampNearZero(a);
            clampNearZero(b);

            wakeIfNeeded(a);
            wakeIfNeeded(b);
        }
    }
}

function resolveBoundaryCollisions(obj) {
    if (obj.isSleeping) return;
    const verts = getPolygonVertices(obj);
    let minX, maxX, minY, maxY;

    if (!verts) {
        minX = obj.x - obj.size;
        maxX = obj.x + obj.size;
        minY = obj.y - obj.size;
        maxY = obj.y + obj.size;
    } else {
        minX = verts[0].x;
        maxX = verts[0].x;
        minY = verts[0].y;
        maxY = verts[0].y;
        for (let i = 1; i < verts.length; i++) {
            const v = verts[i];
            if (v.x < minX) minX = v.x;
            if (v.x > maxX) maxX = v.x;
            if (v.y < minY) minY = v.y;
            if (v.y > maxY) maxY = v.y;
        }
    }

    // Left wall
    if (minX < 0) {
        obj.x += -minX;
        if (obj.vx < 0) obj.vx = -obj.vx * WALL_RESTITUTION;
        obj.angularVelocity += (obj.vx * SURFACE_FRICTION) / Math.max(8, obj.size);
    }

    // Right wall
    if (maxX > canvas.width) {
        obj.x -= (maxX - canvas.width);
        if (obj.vx > 0) obj.vx = -obj.vx * WALL_RESTITUTION;
        obj.angularVelocity -= (obj.vx * SURFACE_FRICTION) / Math.max(8, obj.size);
    }

    // Ceiling
    if (minY < 0) {
        obj.y += -minY;
        if (obj.vy < 0) obj.vy = -obj.vy * FLOOR_RESTITUTION;
    }

    // Floor
    if (maxY > canvas.height) {
        obj.y -= (maxY - canvas.height);
        if (obj.vy > 0) {
            obj.vy = -obj.vy * FLOOR_RESTITUTION;
            if (Math.abs(obj.vy) < STOP_BOUNCE_THRESHOLD) obj.vy = 0;
        }
        obj.vx *= 1 - SURFACE_FRICTION;
        obj.angularVelocity *= 1 - SURFACE_FRICTION;

        // Extra settling to kill micro-jitter when objects are resting on the floor
        if (obj.vy === 0) {
            obj.vx *= FLOOR_SETTLE_DAMPING;
            obj.angularVelocity *= FLOOR_SETTLE_DAMPING;
            if (Math.abs(obj.vx) < 0.03) obj.vx = 0;
            if (Math.abs(obj.angularVelocity) < 0.01) obj.angularVelocity = 0;
        }
    }

    clampNearZero(obj);
}

function clamp(v, min, max) {
    return Math.max(min, Math.min(max, v));
}

function isObjectResting(obj) {
    return (
        Math.abs(obj.vx) < REST_VELOCITY_THRESHOLD &&
        Math.abs(obj.vy) < REST_VELOCITY_THRESHOLD &&
        Math.abs(obj.angularVelocity) < REST_ANGULAR_THRESHOLD
    );
}

function wakeIfNeeded(obj) {
    if (!obj.isSleeping) return;
    const fast =
        Math.abs(obj.vx) > WAKE_VELOCITY_THRESHOLD ||
        Math.abs(obj.vy) > WAKE_VELOCITY_THRESHOLD ||
        Math.abs(obj.angularVelocity) > WAKE_VELOCITY_THRESHOLD;
    if (fast) obj.isSleeping = false;
}

function createHeldObject() {
    const randomType = SPAWN_TYPES[Math.floor(Math.random() * SPAWN_TYPES.length)];
    const obj = new GameObject(0, 0, randomType);
    obj.vx = 0;
    obj.vy = 0;
    obj.angularVelocity = 0;
    obj.angle = 0;
    return obj;
}

function updateHeldObjectPosition() {
    if (!heldObject) return;
    const s = heldObject.size;
    heldObject.x = clamp(heldX, s, canvas.width - s);
    heldObject.y = s + 8;
}

function dropHeldObject() {
    if (!gameRunning || gamePaused || gameOver) return;
    if (!heldObject) return;

    updateHeldObjectPosition();
    heldObject.vx = 0;
    heldObject.vy = 0;
    heldObject.restFrames = 0;
    heldObject.isSleeping = false;
    gameObjects.push(heldObject);
    heldObject = createHeldObject();
    updateHeldObjectPosition();
}

// Object class for falling items
class GameObject {
    constructor(x, y, type, size = 20) {
        this.x = x;
        this.y = y;
        this.type = type;
        this.size = SHAPE_CONFIG[type]?.size ?? size;
        this.vx = (Math.random() - 0.5) * 2; // Random horizontal velocity
        this.vy = 0; // Initial vertical velocity
        this.angle = Math.random() * Math.PI * 2;
        this.angularVelocity = 0;
        this.mass = Math.max(0.5, (this.size * this.size) / 400);
        this.inertia = 0.5 * this.mass * this.size * this.size;
        this.baseColor = OBJECT_TYPES[type].color;
        this.color = this.baseColor;
        this.isColliding = false;
        this.merged = false; // Flag to prevent multiple merges
        this.restFrames = 0;
        this.isSleeping = false;
    }

    // Update object position and velocity
    update() {
        if (this.isSleeping) return;
        // Apply gravity
        this.vy += GRAVITY;

        // Apply friction to horizontal movement
        this.vx *= FRICTION;

        // Update position
        this.x += this.vx;
        this.y += this.vy;

        clampNearZero(this);

        // Disable rotation for stability
        this.angularVelocity = 0;
    }

    // Draw the object on canvas
    draw() {
        ctx.fillStyle = this.color;
        ctx.strokeStyle = '#333';
        ctx.lineWidth = 2;

        switch (this.type) {
            case 'CIRCLE':
                ctx.beginPath();
                ctx.arc(this.x, this.y, this.size, 0, Math.PI * 2);
                ctx.fill();
                ctx.stroke();
                break;

            default: {
                const sides = POLY_SIDES[this.type];
                if (!sides) break;

                ctx.save();
                ctx.translate(this.x, this.y);
                ctx.rotate(this.angle);
                ctx.beginPath();
                for (let i = 0; i < sides; i++) {
                    const a = (Math.PI * 2 * i) / sides - Math.PI / 2;
                    const px = Math.cos(a) * this.size;
                    const py = Math.sin(a) * this.size;
                    if (i === 0) ctx.moveTo(px, py);
                    else ctx.lineTo(px, py);
                }
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
                ctx.restore();
                break;
            }
        }
    }

    // Check collision with another object
    collidesWith(other) {
        const distance = Math.sqrt((this.x - other.x) ** 2 + (this.y - other.y) ** 2);
        return distance < this.size + other.size;
    }

    // Merge with another object of same type
    merge(other) {
        if (this.type === other.type && !this.merged && !other.merged) {
            // Increase size (but cap at maximum)
            this.size = Math.min(this.size + other.size * 0.5, 60);

            // Mark both as merged
            this.merged = true;
            other.merged = true;

            // Add score
            score += Math.floor(this.size);

            // Conservation of momentum
            this.vx = (this.vx + other.vx) / 2;
            this.vy = (this.vy + other.vy) / 2;

            return true;
        }
        return false;
    }
}

// Game loop
function gameLoop() {
    if (!gameRunning || gamePaused) return;

    // Clear canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    if (gameOver) {
        ctx.fillStyle = 'rgba(0, 0, 0, 0.55)';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 48px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('Game Over', canvas.width / 2, canvas.height / 2);

        ctx.font = 'bold 24px Arial';
        ctx.fillText(`Score: ${score}`, canvas.width / 2, canvas.height / 2 + 48);

        ctx.textAlign = 'start';
        ctx.textBaseline = 'alphabetic';
        return;
    }

    // Draw game-over boundary line
    drawTopBoundaryLine();

    // Update all objects
    for (let i = gameObjects.length - 1; i >= 0; i--) {
        const obj = gameObjects[i];

        // Remove merged objects
        if (obj.merged) {
            gameObjects.splice(i, 1);
            continue;
        }

        obj.update();
        resolveBoundaryCollisions(obj);
    }

    // Merge objects that collide and have the same size
    handleSameSizeMerges();

    // Iterative solver: multiple passes to prevent deep interpenetration under pressure
    for (let k = 0; k < SOLVER_ITERATIONS; k++) {
        resolveObjectCollisions();
        for (const obj of gameObjects) {
            resolveBoundaryCollisions(obj);
        }
    }

    // Track rest frames after physics resolution
    for (const obj of gameObjects) {
        if (obj.isSleeping) continue;
        if (isObjectResting(obj)) {
            obj.restFrames += 1;
            if (obj.restFrames >= SLEEP_FRAMES_REQUIRED) {
                obj.vx = 0;
                obj.vy = 0;
                obj.angularVelocity = 0;
                obj.isSleeping = true;
            }
        } else {
            obj.restFrames = 0;
        }
    }

    // Detect circle-vs-circle collisions (no merging yet)
    handleCircleCollisions();

    // Game over if any object crosses the top boundary
    for (const obj of gameObjects) {
        if ((obj.isSleeping || obj.restFrames >= REST_FRAMES_REQUIRED) && obj.y - obj.size < TOP_GAME_OVER_BOUNDARY) {
            gameOver = true;
            break;
        }
    }

    if (gameOver) {
        ctx.fillStyle = 'rgba(0, 0, 0, 0.55)';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 48px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('Game Over', canvas.width / 2, canvas.height / 2);

        ctx.font = 'bold 24px Arial';
        ctx.fillText(`Score: ${score}`, canvas.width / 2, canvas.height / 2 + 48);

        ctx.textAlign = 'start';
        ctx.textBaseline = 'alphabetic';
        return;
    }

    // Draw all objects
    for (let i = 0; i < gameObjects.length; i++) {
        gameObjects[i].draw();
    }

    // Draw held (preview) object last so it's visible
    if (heldObject) {
        updateHeldObjectPosition();
        heldObject.draw();
    }

    // Update score display
    scoreValue.textContent = score;

    // Continue game loop
    animationId = requestAnimationFrame(gameLoop);
}

// Start game
function startGame() {
    if (gameRunning) return;

    gameRunning = true;
    gamePaused = false;
    gameOver = false;
    score = 0;
    gameObjects = [];

    heldX = canvas.width / 2;
    heldObject = createHeldObject();
    updateHeldObjectPosition();

    gameLoop();
}

// Pause game
function pauseGame() {
    gamePaused = !gamePaused;
    pauseBtn.textContent = gamePaused ? 'Resume' : 'Pause';

    if (!gamePaused) {
        gameLoop();
    }
}

// Reset game
function resetGame() {
    gameRunning = false;
    gamePaused = false;
    gameOver = false;
    score = 0;
    gameObjects = [];
    activeCircleCollisions = new Set();
    heldObject = null;
    heldX = canvas.width / 2;

    if (animationId) {
        cancelAnimationFrame(animationId);
    }

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    scoreValue.textContent = score;
    pauseBtn.textContent = 'Pause';
}

// Event listeners
startBtn.addEventListener('click', startGame);
pauseBtn.addEventListener('click', pauseGame);
resetBtn.addEventListener('click', resetGame);

canvas.addEventListener('mousemove', (e) => {
    const rect = canvas.getBoundingClientRect();
    heldX = e.clientX - rect.left;
});

canvas.addEventListener('click', () => {
    dropHeldObject();
});

document.addEventListener('keydown', (e) => {
    if (e.code === 'Space') {
        e.preventDefault();
        dropHeldObject();
    }
});

// Initialize game
resetGame();
