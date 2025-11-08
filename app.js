// Utilidades generales -------------------------------------------------------
const tabs = document.querySelectorAll(".tab-button");
const sections = document.querySelectorAll("[role=tabpanel]");

const randomExamples = [
    "x^2 + y^2",
    "sin(x) * cos(y)",
    "exp(-x^2 - y^2)",
    "x * y + z",
    "sin(x) + y^3 - 2"
];

const formatNumber = (value, digits = 4) => {
    if (!Number.isFinite(value)) return "NaN";
    const rounded = Number(value.toFixed(digits));
    return Math.abs(rounded) < 1e-9 ? "0" : rounded.toString();
};

const linspace = (min, max, count) => {
    if (count <= 1) return [min];
    const step = (max - min) / (count - 1);
    return Array.from({ length: count }, (_, i) => min + step * i);
};

const detectVariables = (expression) => {
    try {
        const node = math.parse(expression);
        const vars = new Set();
        node.traverse((child) => {
            if (child.isSymbolNode && ["x", "y", "z"].includes(child.name)) {
                vars.add(child.name);
            }
        });
        return Array.from(vars).sort();
    } catch (error) {
        return ["x", "y"];
    }
};

const evaluateCompiled = (compiled, scope) => {
    try {
        const value = compiled.evaluate(scope);
        return Number.isFinite(value) ? value : null;
    } catch (error) {
        return null;
    }
};

const analyzeGrid = (grid, xValues, yValues) => {
    let minZ = Infinity;
    let maxZ = -Infinity;
    let minX = Infinity;
    let maxX = -Infinity;
    let minY = Infinity;
    let maxY = -Infinity;
    let validCount = 0;

    grid.forEach((row, j) => {
        row.forEach((value, i) => {
            if (value !== null) {
                const x = xValues[i];
                const y = yValues[j];
                minX = Math.min(minX, x);
                maxX = Math.max(maxX, x);
                minY = Math.min(minY, y);
                maxY = Math.max(maxY, y);
                minZ = Math.min(minZ, value);
                maxZ = Math.max(maxZ, value);
                validCount += 1;
            }
        });
    });

    return {
        hasData: validCount > 0,
        count: validCount,
        x: [minX, maxX],
        y: [minY, maxY],
        z: [minZ, maxZ]
    };
};

const showDomainRange = (stats, zParam = null) => {
    const domainNode = document.getElementById("domain-output");
    const rangeNode = document.getElementById("range-output");

    if (!stats.hasData) {
        domainNode.textContent = "Sin puntos válidos en el rango evaluado.";
        rangeNode.textContent = "-";
        return;
    }

    domainNode.textContent = `x ? [${formatNumber(stats.x[0])}, ${formatNumber(stats.x[1])}], y ? [${formatNumber(stats.y[0])}, ${formatNumber(stats.y[1])}]${zParam !== null ? `, z fijo = ${formatNumber(zParam)}` : ""}`;
    rangeNode.textContent = `f(x, y) ? [${formatNumber(stats.z[0])}, ${formatNumber(stats.z[1])}]`;
};

// Tabs -----------------------------------------------------------------------
tabs.forEach((button) => {
    button.addEventListener("click", () => {
        const target = button.dataset.target;
        tabs.forEach((btn) => btn.setAttribute("aria-selected", btn === button ? "true" : "false"));
        sections.forEach((section) => {
            section.classList.toggle("hidden", section.id !== target);
        });
    });
});

// Visualización --------------------------------------------------------------
const plotButton = document.getElementById("plot-btn");
const randomButton = document.getElementById("random-example-btn");

const plotFunction = () => {
    const expression = document.getElementById("function-input").value.trim();
    const resolution = Math.max(10, Number(document.getElementById("resolution").value) || 35);
    const xMin = Number(document.getElementById("x-min").value) || -5;
    const xMax = Number(document.getElementById("x-max").value) || 5;
    const yMin = Number(document.getElementById("y-min").value) || -5;
    const yMax = Number(document.getElementById("y-max").value) || 5;
    const zParam = Number(document.getElementById("z-param").value) || 0;

    try {
        const compiled = math.compile(expression);
        const variables = detectVariables(expression);
        const hasZ = variables.includes("z");
        const xValues = linspace(xMin, xMax, resolution);
        const yValues = linspace(yMin, yMax, resolution);

        const grid = yValues.map((y) => {
            return xValues.map((x) => {
                const scope = { x, y };
                if (hasZ) scope.z = zParam;
                return evaluateCompiled(compiled, scope);
            });
        });

        const stats = analyzeGrid(grid, xValues, yValues);
        showDomainRange(stats, hasZ ? zParam : null);

        const data = [{
            type: "surface",
            x: xValues,
            y: yValues,
            z: grid,
            colorscale: "Viridis",
            showscale: true,
            hovertemplate: "x=%{x:.2f}<br>y=%{y:.2f}<br>f=%{z:.2f}<extra></extra>"
        }];

        Plotly.newPlot("plot-visual", data, {
            title: hasZ ? `Gráfica de f(x, y, z=${formatNumber(zParam, 2)})` : "Gráfica de f(x, y)",
            paper_bgcolor: "rgba(0,0,0,0)",
            plot_bgcolor: "rgba(0,0,0,0)",
            scene: {
                xaxis: { title: "x" },
                yaxis: { title: "y" },
                zaxis: { title: "f(x, y)" }
            }
        }, { responsive: true });
    } catch (error) {
        alert("No se pudo graficar la función. Revisa la sintaxis.");
        console.error(error);
    }
};

plotButton.addEventListener("click", plotFunction);
randomButton.addEventListener("click", () => {
    const choice = randomExamples[Math.floor(Math.random() * randomExamples.length)];
    document.getElementById("function-input").value = choice;
    document.getElementById("derivative-function").value = choice;
    plotFunction();
});

// Derivadas y gradiente ------------------------------------------------------
const derivativeButton = document.getElementById("derivative-btn");

const showDerivatives = () => {
    const expression = document.getElementById("derivative-function").value.trim();
    const x0 = Number(document.getElementById("point-x").value) || 0;
    const y0 = Number(document.getElementById("point-y").value) || 0;
    const z0 = Number(document.getElementById("point-z").value) || 0;

    try {
        const variables = detectVariables(expression);
        const compiled = math.compile(expression);
        const partials = variables.map((variable) => ({
            variable,
            derivative: math.derivative(expression, variable).toString()
        }));

        const partialsText = partials.map((p) => `?f/?${p.variable} = ${p.derivative}`).join("\n");
        document.getElementById("partials-output").textContent = partialsText || "-";

        const scope = { x: x0, y: y0, z: z0 };
        const gradientValues = partials.map((p) => {
            const compiledPartial = math.compile(p.derivative);
            return evaluateCompiled(compiledPartial, scope);
        });

        const gradientText = `?f(${x0}, ${y0}${variables.includes("z") ? ", " + z0 : ""}) = ?${gradientValues.map((val) => formatNumber(val || 0)).join(", ")}?`;
        document.getElementById("gradient-output").textContent = gradientText;
    } catch (error) {
        document.getElementById("partials-output").textContent = "No se pudo derivar la función.";
        document.getElementById("gradient-output").textContent = "-";
        console.error(error);
    }
};

derivativeButton.addEventListener("click", showDerivatives);

// Integrales múltiples -------------------------------------------------------
const integralButton = document.getElementById("integral-btn");

const evaluateIntegrand = (compiled, scope) => {
    const value = evaluateCompiled(compiled, scope);
    return value ?? 0;
};

const approximateDoubleIntegral = (compiled, bounds, subdivisions) => {
    const { x, y } = bounds;
    const dx = (x[1] - x[0]) / subdivisions;
    const dy = (y[1] - y[0]) / subdivisions;
    let sum = 0;
    for (let i = 0; i < subdivisions; i++) {
        for (let j = 0; j < subdivisions; j++) {
            const xi = x[0] + dx * (i + 0.5);
            const yj = y[0] + dy * (j + 0.5);
            sum += evaluateIntegrand(compiled, { x: xi, y: yj });
        }
    }
    return sum * dx * dy;
};

const approximateTripleIntegral = (compiled, bounds, subdivisions) => {
    const { x, y, z } = bounds;
    const dx = (x[1] - x[0]) / subdivisions;
    const dy = (y[1] - y[0]) / subdivisions;
    const dz = (z[1] - z[0]) / subdivisions;
    let sum = 0;
    for (let i = 0; i < subdivisions; i++) {
        for (let j = 0; j < subdivisions; j++) {
            for (let k = 0; k < subdivisions; k++) {
                const xi = x[0] + dx * (i + 0.5);
                const yj = y[0] + dy * (j + 0.5);
                const zk = z[0] + dz * (k + 0.5);
                sum += evaluateIntegrand(compiled, { x: xi, y: yj, z: zk });
            }
        }
    }
    return sum * dx * dy * dz;
};

const renderIntegralPlot = (type, compiled, bounds) => {
    try {
        if (type === "double") {
            const xVals = linspace(bounds.x[0], bounds.x[1], 30);
            const yVals = linspace(bounds.y[0], bounds.y[1], 30);
            const grid = yVals.map((y) => xVals.map((x) => evaluateCompiled(compiled, { x, y })));
            Plotly.newPlot("plot-integral", [{
                type: "surface",
                x: xVals,
                y: yVals,
                z: grid,
                colorscale: "Magma",
                showscale: true
            }], {
                title: "Integrando sobre la región",
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)"
            }, { responsive: true });
        } else {
            const samples = 25;
            const points = [];
            for (let i = 0; i < samples; i++) {
                const sx = bounds.x[0] + Math.random() * (bounds.x[1] - bounds.x[0]);
                const sy = bounds.y[0] + Math.random() * (bounds.y[1] - bounds.y[0]);
                const sz = bounds.z[0] + Math.random() * (bounds.z[1] - bounds.z[0]);
                points.push({ x: sx, y: sy, z: sz, value: evaluateCompiled(compiled, { x: sx, y: sy, z: sz }) });
            }
            Plotly.newPlot("plot-integral", [{
                type: "scatter3d",
                mode: "markers",
                x: points.map((p) => p.x),
                y: points.map((p) => p.y),
                z: points.map((p) => p.z),
                marker: {
                    size: 5,
                    color: points.map((p) => p.value),
                    colorscale: "Viridis",
                    colorbar: { title: "f(x,y,z)" }
                }
            }], {
                title: "Muestreo del integrando",
                paper_bgcolor: "rgba(0,0,0,0)",
                plot_bgcolor: "rgba(0,0,0,0)"
            }, { responsive: true });
        }
    } catch (error) {
        console.error("No se pudo graficar el integrando", error);
    }
};

const approximateIntegral = () => {
    const expression = document.getElementById("integral-function").value.trim();
    const type = document.getElementById("integral-type").value;
    const subdivisions = Math.max(5, Number(document.getElementById("integral-resolution").value) || 20);
    const bounds = {
        x: [Number(document.getElementById("ix-min").value) || 0, Number(document.getElementById("ix-max").value) || 1],
        y: [Number(document.getElementById("iy-min").value) || 0, Number(document.getElementById("iy-max").value) || 1],
        z: [Number(document.getElementById("iz-min").value) || 0, Number(document.getElementById("iz-max").value) || 1]
    };

    try {
        const compiled = math.compile(expression);
        const result = type === "double"
            ? approximateDoubleIntegral(compiled, bounds, subdivisions)
            : approximateTripleIntegral(compiled, bounds, subdivisions);

        document.getElementById("integral-output").textContent = `˜ ${formatNumber(result, 6)}`;
        renderIntegralPlot(type, compiled, bounds);
    } catch (error) {
        document.getElementById("integral-output").textContent = "No se pudo evaluar la integral.";
        console.error(error);
    }
};

integralButton.addEventListener("click", approximateIntegral);

// Optimización con restricciones --------------------------------------------
const optimizeButton = document.getElementById("optimize-btn");

const gaussianSolve = (matrix, vector) => {
    const n = matrix.length;
    const augmented = matrix.map((row, i) => [...row, vector[i]]);

    for (let i = 0; i < n; i++) {
        let pivot = i;
        for (let j = i + 1; j < n; j++) {
            if (Math.abs(augmented[j][i]) > Math.abs(augmented[pivot][i])) {
                pivot = j;
            }
        }
        if (Math.abs(augmented[pivot][i]) < 1e-9) throw new Error("Jacobiano singular");
        [augmented[i], augmented[pivot]] = [augmented[pivot], augmented[i]];

        const pivotValue = augmented[i][i];
        for (let k = i; k <= n; k++) {
            augmented[i][k] /= pivotValue;
        }

        for (let r = 0; r < n; r++) {
            if (r !== i) {
                const factor = augmented[r][i];
                for (let c = i; c <= n; c++) {
                    augmented[r][c] -= factor * augmented[i][c];
                }
            }
        }
    }

    return augmented.map((row) => row[n]);
};

const numericalJacobian = (fn, vars) => {
    const h = 1e-4;
    const base = fn(vars);
    return base.map(() => Array(vars.length).fill(0)).map((row, rowIdx) => {
        return row.map((_, colIdx) => {
            const forward = [...vars];
            forward[colIdx] += h;
            const next = fn(forward);
            return (next[rowIdx] - base[rowIdx]) / h;
        });
    });
};

const newtonSystem = (fn, initialGuess, maxIterations = 30, tolerance = 1e-6) => {
    let vars = [...initialGuess];
    for (let iteration = 0; iteration < maxIterations; iteration++) {
        const F = fn(vars);
        const norm = Math.sqrt(F.reduce((acc, value) => acc + value * value, 0));
        if (norm < tolerance) return vars;
        const J = numericalJacobian(fn, vars);
        const delta = gaussianSolve(J, F.map((value) => -value));
        vars = vars.map((value, idx) => value + delta[idx]);
    }
    throw new Error("Newton no convergió");
};

const solveLagrange = ({ objective, constraint, constant, initial }) => {
    const fx = math.derivative(objective, "x").compile();
    const fy = math.derivative(objective, "y").compile();
    const gx = math.derivative(constraint, "x").compile();
    const gy = math.derivative(constraint, "y").compile();
    const g = math.compile(constraint);
    const f = math.compile(objective);

    const system = (vars) => {
        const [x, y, lambda] = vars;
        const scope = { x, y };
        return [
            fx.evaluate(scope) - lambda * gx.evaluate(scope),
            fy.evaluate(scope) - lambda * gy.evaluate(scope),
            g.evaluate(scope) - constant
        ];
    };

    const solution = newtonSystem(system, [initial.x, initial.y, 0]);
    const [x, y, lambda] = solution;
    return { x, y, lambda, value: f.evaluate({ x, y }) };
};

const renderOptimizationPlot = (objective, constraint, constant, solution) => {
    try {
        const compiledF = math.compile(objective);
        const compiledG = math.compile(constraint);
        const span = 3;
        const xRange = [solution ? solution.x - span : -2, solution ? solution.x + span : 2];
        const yRange = [solution ? solution.y - span : -2, solution ? solution.y + span : 2];
        const xVals = linspace(xRange[0], xRange[1], 60);
        const yVals = linspace(yRange[0], yRange[1], 60);

        const surface = yVals.map((y) => xVals.map((x) => evaluateCompiled(compiledF, { x, y })));
        const constraintGrid = yVals.map((y) => xVals.map((x) => evaluateCompiled(compiledG, { x, y }) - constant));

        const data = [
            {
                type: "heatmap",
                x: xVals,
                y: yVals,
                z: surface,
                colorscale: "Electric",
                name: "f(x, y)",
                showscale: true
            },
            {
                type: "contour",
                x: xVals,
                y: yVals,
                z: constraintGrid,
                contours: {
                    start: 0,
                    end: 0,
                    size: 0.5,
                    coloring: "lines"
                },
                line: { color: "#f97316", width: 3 },
                name: "g(x, y) = c"
            }
        ];

        if (solution) {
            data.push({
                type: "scatter",
                x: [solution.x],
                y: [solution.y],
                mode: "markers",
                marker: { size: 12, color: "#22d3ee" },
                name: "Óptimo"
            });
        }

        Plotly.newPlot("plot-optimization", data, {
            title: "Nivel de f(x, y) y restricción",
            paper_bgcolor: "rgba(0,0,0,0)",
            plot_bgcolor: "rgba(0,0,0,0)",
            xaxis: { title: "x" },
            yaxis: { title: "y" }
        }, { responsive: true });
    } catch (error) {
        console.error("No se pudo graficar la optimización", error);
    }
};

const runOptimization = () => {
    const objective = document.getElementById("objective-input").value.trim();
    const constraint = document.getElementById("constraint-input").value.trim();
    const constant = Number(document.getElementById("constraint-value").value) || 0;
    const initial = {
        x: Number(document.getElementById("opt-x0").value) || 0,
        y: Number(document.getElementById("opt-y0").value) || 0
    };

    try {
        const solution = solveLagrange({ objective, constraint, constant, initial });
        document.getElementById("optimum-output").textContent = `(${formatNumber(solution.x)}, ${formatNumber(solution.y)})`;
        document.getElementById("lambda-output").textContent = formatNumber(solution.lambda);
        document.getElementById("optimum-value-output").textContent = formatNumber(solution.value);
        renderOptimizationPlot(objective, constraint, constant, solution);
    } catch (error) {
        document.getElementById("optimum-output").textContent = "No convergió";
        document.getElementById("lambda-output").textContent = "-";
        document.getElementById("optimum-value-output").textContent = "-";
        alert("El método de Newton no convergió. Intenta con otros valores iniciales.");
        console.error(error);
        renderOptimizationPlot(objective, constraint, constant, null);
    }
};

optimizeButton.addEventListener("click", runOptimization);

// Inicialización -------------------------------------------------------------
plotFunction();
showDerivatives();
renderOptimizationPlot(
    document.getElementById("objective-input").value,
    document.getElementById("constraint-input").value,
    Number(document.getElementById("constraint-value").value) || 0,
    null
);
