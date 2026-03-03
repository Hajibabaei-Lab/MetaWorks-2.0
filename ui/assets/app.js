const API_BASE = window.location.origin;

const STORAGE_KEYS = {
  trackedRuns: "metaworks.trackedRuns",
  draft: "metaworks.draft",
  overridesPrefix: "metaworks.overrides.",
};

const BUNDLED_TEST = {
  runtime: "docker",
  containerImage: "metaworks:latest",
  inputDir: "/workspace/tests/testing_data",
  adapterPath: "/workspace/tests/adapters_anchored.fasta",
};

function nowIso() {
  return new Date().toISOString();
}

function safeJsonParse(text) {
  try {
    return { ok: true, value: JSON.parse(text) };
  } catch (err) {
    return { ok: false, error: err };
  }
}

function storageGet(key, fallback) {
  try {
    const raw = localStorage.getItem(key);
    if (!raw) return fallback;
    return JSON.parse(raw);
  } catch {
    return fallback;
  }
}

function storageSet(key, value) {
  localStorage.setItem(key, JSON.stringify(value));
}

function el(tag, attrs = {}, children = []) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "class") node.className = v;
    else if (k === "text") node.textContent = v;
    else if (k === "html") node.innerHTML = v;
    else if (k.startsWith("on") && typeof v === "function")
      node.addEventListener(k.slice(2).toLowerCase(), v);
    else if (v === true) node.setAttribute(k, "");
    else if (v !== false && v != null) node.setAttribute(k, String(v));
  }
  for (const child of Array.isArray(children) ? children : [children]) {
    if (child == null) continue;
    node.appendChild(typeof child === "string" ? document.createTextNode(child) : child);
  }
  return node;
}

function toast(title, message) {
  const root = document.getElementById("toasts");
  const node = el("div", { class: "toast", role: "status" }, [
    el("div", { class: "toast-title", text: title }),
    el("div", { class: "toast-body", text: message }),
  ]);
  root.appendChild(node);
  setTimeout(() => node.remove(), 5200);
}

async function apiFetch(path, options = {}) {
  const res = await fetch(`${API_BASE}${path}`, options);
  const contentType = res.headers.get("content-type") || "";
  const isJson = contentType.includes("application/json");
  const body = isJson ? await res.json().catch(() => null) : await res.text().catch(() => "");
  if (!res.ok) {
    const detail =
      (body && typeof body === "object" && body.detail) ||
      (typeof body === "string" && body) ||
      `HTTP ${res.status}`;
    throw new Error(detail);
  }
  return body;
}

const api = {
  health: () => apiFetch("/health", { headers: { Accept: "text/plain" } }),
  getSettingsPaths: () => apiFetch("/settings/paths"),
  listProfiles: () => apiFetch("/profiles"),
  getProfile: (name) => apiFetch(`/profiles/${encodeURIComponent(name)}`),
  submitRun: (payload) =>
    apiFetch("/runs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    }),
  getRunStatus: (runId) => apiFetch(`/runs/${encodeURIComponent(runId)}`),
  getRunLogs: (runId, lines = 200) =>
    apiFetch(`/runs/${encodeURIComponent(runId)}/logs?lines=${encodeURIComponent(lines)}`),
  cancelRun: (runId) =>
    apiFetch(`/runs/${encodeURIComponent(runId)}/cancel`, { method: "POST" }),
  deleteRun: (runId) => apiFetch(`/runs/${encodeURIComponent(runId)}`, { method: "DELETE" }),
  getConfigDefaultsText: (workflow) => apiFetch(`/configs/defaults/${workflow}`),
  getConfigDefaultsSections: (workflow) =>
    apiFetch(`/configs/defaults/${workflow}/sections`),
  renderConfig: (profile, workflow, overrides) =>
    apiFetch(`/configs/render`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ profile, workflow, overrides }),
    }),
  listAdapters: () => apiFetch("/adapters"),
  listClassifiers: () => apiFetch("/classifiers"),
  deleteAdapter: (name) =>
    apiFetch(`/delete/adapters/${encodeURIComponent(name)}`, { method: "POST" }),
  deleteClassifier: (name) =>
    apiFetch(`/delete/classifiers/${encodeURIComponent(name)}`, { method: "POST" }),
};

function deepEqual(a, b) {
  if (a === b) return true;
  if (typeof a !== typeof b) return false;
  if (a == null || b == null) return false;
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!deepEqual(a[i], b[i])) return false;
    return true;
  }
  if (typeof a === "object") {
    const ak = Object.keys(a);
    const bk = Object.keys(b);
    if (ak.length !== bk.length) return false;
    for (const k of ak) if (!deepEqual(a[k], b[k])) return false;
    return true;
  }
  return false;
}

function getAt(obj, path) {
  let cur = obj;
  for (const key of path) {
    if (cur == null) return undefined;
    cur = cur[key];
  }
  return cur;
}

function setAt(obj, path, value) {
  let cur = obj;
  for (let i = 0; i < path.length - 1; i++) {
    const key = path[i];
    if (!cur[key] || typeof cur[key] !== "object") cur[key] = {};
    cur = cur[key];
  }
  cur[path[path.length - 1]] = value;
}

function deleteAt(obj, path) {
  const stack = [];
  let cur = obj;
  for (let i = 0; i < path.length - 1; i++) {
    const key = path[i];
    if (!cur[key] || typeof cur[key] !== "object") return;
    stack.push([cur, key]);
    cur = cur[key];
  }
  delete cur[path[path.length - 1]];
  for (let i = stack.length - 1; i >= 0; i--) {
    const [parent, key] = stack[i];
    if (parent[key] && typeof parent[key] === "object" && Object.keys(parent[key]).length === 0) {
      delete parent[key];
    }
  }
}

function mergeViewValue(defaultValue, overrideValue) {
  if (overrideValue === undefined) return defaultValue;
  return overrideValue;
}

function statusPill(status) {
  const normalized = String(status || "").toLowerCase();
  const cls =
    normalized === "completed"
      ? "ok"
      : normalized === "failed"
        ? "bad"
        : normalized === "cancelled"
          ? "warn"
          : normalized === "running"
            ? "neutral"
            : normalized === "queued"
              ? "neutral"
              : normalized === "staged"
                ? "neutral"
                : "neutral";
  const label = normalized || "unknown";
  return el("span", { class: `pill ${cls}`, text: label });
}

function progressForStatus(status) {
  const s = String(status || "").toLowerCase();
  if (s === "queued") return { kind: "determinate", pct: 8 };
  if (s === "staged") return { kind: "determinate", pct: 12 };
  if (s === "running") return { kind: "indeterminate", pct: 0 };
  if (s === "completed") return { kind: "determinate", pct: 100 };
  if (s === "failed") return { kind: "determinate", pct: 100 };
  if (s === "cancelled") return { kind: "determinate", pct: 100 };
  return { kind: "determinate", pct: 0 };
}

function formatTime(ts) {
  if (!ts) return "—";
  const d = new Date(ts);
  if (Number.isNaN(d.getTime())) return String(ts);
  return d.toLocaleString();
}

const state = {
  draft: storageGet(STORAGE_KEYS.draft, {
    profile: "coi",
    workflow: "esv",
    run_name: "",
    runtime: BUNDLED_TEST.runtime,
    container_image: BUNDLED_TEST.containerImage,
    bind_paths_text: "",
    cache_dir: "",
    input_dir: BUNDLED_TEST.inputDir,
    sample_source: "folder",
    samples_csv: "",
    cores: "",
    mem_gb: "",
    dry_run: false,
    keep_outputs: true,
    notes: "",
  }),
  profiles: [],  // Loaded from API
  overridesByWorkflow: {
    esv: storageGet(`${STORAGE_KEYS.overridesPrefix}esv`, {}),
    otu: storageGet(`${STORAGE_KEYS.overridesPrefix}otu`, {}),
  },
  defaultsByWorkflow: { esv: null, otu: null },
  systemSettings: {
    allowed_runtimes: ["docker", "apptainer"],
    retention_policy: "until_download",
  },
  trackedRuns: storageGet(STORAGE_KEYS.trackedRuns, []),
  runStatusById: {},
  view: { route: "submit", intervalId: null, runsPollInFlight: false },
  lastHealth: null,
};

function normalizeDraftRuntime() {
  const allowed = Array.isArray(state.systemSettings.allowed_runtimes)
    ? state.systemSettings.allowed_runtimes
    : ["docker", "apptainer"];
  const fallback = allowed[0] || "docker";
  if (!allowed.includes(state.draft.runtime)) {
    state.draft.runtime = fallback;
    persistDraft();
  }
}

function isActiveRunStatus(status) {
  return ["running", "queued", "staged"].includes(String(status || "").toLowerCase());
}

async function pollActiveRunsAndMaybeRender() {
  if (state.view.route !== "runs") return;
  if (state.view.runsPollInFlight) return;
  if (!state.trackedRuns.length) return;

  const activeRunIds = state.trackedRuns.filter((runId) => {
    const known = state.runStatusById?.[runId];
    return !known || isActiveRunStatus(known.status);
  });

  if (!activeRunIds.length) return;

  state.view.runsPollInFlight = true;
  try {
    const updates = await Promise.all(
      activeRunIds.map(async (runId) => {
        try {
          const next = await api.getRunStatus(runId);
          return { runId, next };
        } catch {
          return null;
        }
      })
    );

    let changed = false;
    state.runStatusById = state.runStatusById || {};

    for (const row of updates) {
      if (!row) continue;
      const prev = state.runStatusById[row.runId];
      state.runStatusById[row.runId] = row.next;
      if (!prev || !deepEqual(prev, row.next)) changed = true;
    }

    if (changed && state.view.route === "runs") render();
  } finally {
    state.view.runsPollInFlight = false;
  }
}

function persistDraft() {
  storageSet(STORAGE_KEYS.draft, state.draft);
}

function persistOverrides(workflow) {
  storageSet(`${STORAGE_KEYS.overridesPrefix}${workflow}`, state.overridesByWorkflow[workflow]);
}

function persistTrackedRuns() {
  const unique = Array.from(new Set(state.trackedRuns)).filter(Boolean);
  state.trackedRuns = unique;
  storageSet(STORAGE_KEYS.trackedRuns, unique);
}

function currentOverrides() {
  return state.overridesByWorkflow[state.draft.workflow] || {};
}

function syncOverridesTextarea() {
  const textArea = document.getElementById("overrides_json");
  if (!textArea) return;
  const workflow = state.draft.workflow;
  textArea.value = JSON.stringify(state.overridesByWorkflow[workflow] || {}, null, 2);
}

const MODULE_DEFS = [
  {
    key: "preprocessing",
    label: "Preprocessing",
    help: "Quality filtering and read merging (SeqPrep).",
  },
  {
    key: "trimming",
    label: "Trimming",
    help: "Adapter trimming (Cutadapt).",
  },
  {
    key: "denoising",
    label: "Denoising",
    help: "ESV/OTU inference (VSEARCH).",
  },
  {
    key: "classification",
    label: "Classification",
    help: "Taxonomic assignment (RDP or SINTAX).",
  },
  {
    key: "pseudogene_filtering",
    label: "Pseudogene filtering",
    help: "Optional pseudogene removal (HMM/ORF).",
  },
  {
    key: "stats",
    label: "Stats",
    help: "Summary statistics and reports.",
  },
];

function moduleDefaultsFor(workflow) {
  const fromApi = state.defaultsByWorkflow[workflow];
  if (fromApi && fromApi.modules && typeof fromApi.modules === "object") return fromApi.modules;
  return {
    preprocessing: true,
    trimming: true,
    denoising: true,
    classification: true,
    pseudogene_filtering: false,
    stats: true,
  };
}

function effectiveModuleValue(workflow, moduleKey) {
  const base = moduleDefaultsFor(workflow);
  const overrides = state.overridesByWorkflow[workflow] || {};
  const overrideValue = getAt(overrides, ["modules", moduleKey]);
  return overrideValue === undefined ? Boolean(base[moduleKey]) : Boolean(overrideValue);
}

function setModuleOverride(workflow, moduleKey, enabled) {
  const overrides = state.overridesByWorkflow[workflow] || {};
  const base = moduleDefaultsFor(workflow);
  const defaultEnabled = Boolean(base[moduleKey]);
  const nextEnabled = Boolean(enabled);
  if (nextEnabled === defaultEnabled) deleteAt(overrides, ["modules", moduleKey]);
  else setAt(overrides, ["modules", moduleKey], nextEnabled);
  state.overridesByWorkflow[workflow] = overrides;
  persistOverrides(workflow);
}

function applyPreset(preset, workflow) {
  const w = workflow || state.draft.workflow;
  if (preset === "custom") return {};
  if (preset === "16s") {
    return {
      pipeline: { name: w, output_dir: "16S_results" },
      classification: { marker: "16S" },
    };
  }
  return {
    pipeline: { name: w, output_dir: "COI_results" },
    classification: {
      marker: "COI",
      rdp: { use_custom_classifier: true, classifier_path: "runtime/classifiers/COI.properties" },
    },
  };
}

async function ensureDefaults(workflow) {
  if (state.defaultsByWorkflow[workflow]) return state.defaultsByWorkflow[workflow];
  const res = await api.getConfigDefaultsSections(workflow);
  state.defaultsByWorkflow[workflow] = res.sections || {};
  return state.defaultsByWorkflow[workflow];
}

function renderPageHeader(title, subtitle, right = null) {
  return el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h2", { class: "card-title", text: title }),
        el("div", { class: "muted", text: subtitle }),
      ]),
      right || el("div"),
    ]),
  ]);
}

function renderSubmitView() {
  const root = el("div", { class: "grid", role: "region", "aria-label": "Submit run view" });

  const responseBox = el("pre", { class: "code", text: "" });
  const responseCard = el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [el("h3", { class: "card-title", text: "Response" })]),
      el("div", { class: "row" }, [
        el(
          "button",
          {
            class: "btn btn-secondary",
            type: "button",
            onclick: async () => {
              try {
                const ok = await api.health();
                state.lastHealth = { ok, at: nowIso() };
                toast("Backend", "Health check ok");
              } catch (err) {
                toast("Backend", `Health check failed: ${String(err.message || err)}`);
              }
            },
          },
          "Health"
        ),
      ]),
    ]),
    el("div", { class: "card-body" }, [
      el("div", { class: "help" }, [
        "Successful submissions return a ",
        el("span", { class: "kbd", text: "run_id" }),
        " you can track in the Runs tab.",
      ]),
      el("div", { class: "hr" }),
      responseBox,
    ]),
  ]);

  const form = el("form", {
    onsubmit: async (e) => {
      e.preventDefault();
      responseBox.textContent = "";

      const workflow = state.draft.workflow;
      const runtime = state.draft.runtime;
      const runName = String(state.draft.run_name || "").trim();
      const inputDir = String(state.draft.input_dir || "").trim();
      const sampleSource = state.draft.sample_source;
      const samplesCsv = String(state.draft.samples_csv || "").trim();

      if (!runName) {
        toast("Validation", "Run name is required.");
        return;
      }
      if (!inputDir) {
        toast("Validation", "Input directory is required.");
        return;
      }
      if (sampleSource === "csv" && !samplesCsv) {
        toast("Validation", "Samples CSV is required when sample source is CSV.");
        return;
      }
      if ((runtime === "docker" || runtime === "apptainer") && !state.draft.container_image.trim()) {
        toast("Validation", "Container image is required for Docker/Apptainer runtime.");
        return;
      }

      const bindPaths = String(state.draft.bind_paths_text || "")
        .split("\n")
        .map((s) => s.trim())
        .filter(Boolean);

      const overrides = currentOverrides();

      const payload = {
        profile: state.draft.profile || "coi",
        workflow,
        run_name: runName,
        runtime,
        container_image:
          runtime === "conda" ? null : String(state.draft.container_image || "").trim() || null,
        bind_paths: bindPaths,
        cache_dir: String(state.draft.cache_dir || "").trim() || null,
        config_overrides: overrides,
        input_dir: inputDir,
        sample_source: sampleSource,
        samples_csv: sampleSource === "csv" ? samplesCsv : null,
        notes: String(state.draft.notes || "").trim() || null,
        cores: state.draft.cores ? Number(state.draft.cores) : null,
        mem_gb: state.draft.mem_gb ? Number(state.draft.mem_gb) : null,
        dry_run: Boolean(state.draft.dry_run),
        keep_outputs: true,
      };

      const btn = form.querySelector('button[type="submit"]');
      btn.disabled = true;
      btn.textContent = "Submitting…";
      try {
        const res = await api.submitRun(payload);
        responseBox.textContent = JSON.stringify(res, null, 2);
        if (res && res.run_id) {
          state.trackedRuns.unshift(res.run_id);
          persistTrackedRuns();
          toast("Run submitted", `Tracking ${res.run_id}`);
        } else {
          toast("Submitted", "Run submitted (no run_id returned).");
        }
      } catch (err) {
        responseBox.textContent = `Error: ${String(err.message || err)}`;
        toast("Submit failed", String(err.message || err));
      } finally {
        btn.disabled = false;
        btn.textContent = "Submit Run";
      }
    },
  });

  function rerender() {
    render();
  }

  function applyBundledTestRunDefaults() {
    const ts = new Date();
    const stamp = [
      ts.getFullYear(),
      String(ts.getMonth() + 1).padStart(2, "0"),
      String(ts.getDate()).padStart(2, "0"),
      String(ts.getHours()).padStart(2, "0"),
      String(ts.getMinutes()).padStart(2, "0"),
      String(ts.getSeconds()).padStart(2, "0"),
    ].join("");
    const workflow = state.draft.workflow || "esv";
    state.draft.runtime = BUNDLED_TEST.runtime;
    state.draft.container_image = BUNDLED_TEST.containerImage;
    state.draft.input_dir = BUNDLED_TEST.inputDir;
    state.draft.sample_source = "folder";
    state.draft.samples_csv = "";
    state.draft.keep_outputs = true;
    if (!String(state.draft.run_name || "").trim()) {
      state.draft.run_name = `${workflow}-test-${stamp}`;
    }
    const overrides = state.overridesByWorkflow[workflow] || {};
    setAt(overrides, ["trimming", "adapters"], BUNDLED_TEST.adapterPath);
    setAt(overrides, ["classification", "rdp", "use_custom_classifier"], false);
    setAt(overrides, ["classification", "rdp", "builtin_classifier"], "fungallsu");
    state.overridesByWorkflow[workflow] = overrides;
    persistOverrides(workflow);
    persistDraft();
    toast("Submit", "Bundled Docker test run defaults applied.");
    rerender();
  }

  function inputField({ label, help, input }) {
    return el("div", { class: "field" }, [
      el("label", { for: input.id, text: label }),
      input,
      help ? el("div", { class: "help", text: help }) : null,
    ]);
  }

  const allowedRuntimes = Array.isArray(state.systemSettings.allowed_runtimes)
    ? state.systemSettings.allowed_runtimes
    : ["docker", "apptainer"];
  const retentionPolicy = String(state.systemSettings.retention_policy || "until_download");

  const presetSelect = el(
    "select",
    {
      id: "preset",
      onchange: (e) => {
        const preset = e.target.value;
        state.draft.preset = preset;
        const w = state.draft.workflow;
        state.overridesByWorkflow[w] = applyPreset(preset, w);
        persistOverrides(w);
        persistDraft();
        rerender();
      },
    },
    [
      el("option", { value: "coi", text: "COI Standard" }),
      el("option", { value: "16s", text: "16S Microbiome" }),
      el("option", { value: "custom", text: "Custom (blank)" }),
    ]
  );
  presetSelect.value = state.draft.preset || "coi";

  const workflowSelect = el(
    "select",
    {
      id: "workflow",
      onchange: (e) => {
        const w = e.target.value;
        state.draft.workflow = w;
        if (!state.overridesByWorkflow[w] || Object.keys(state.overridesByWorkflow[w]).length === 0) {
          state.overridesByWorkflow[w] = applyPreset(state.draft.preset, w);
          persistOverrides(w);
        }
        persistDraft();
        rerender();
      },
    },
    [
      el("option", { value: "esv", text: "ESV" }),
      el("option", { value: "otu", text: "OTU" }),
    ]
  );
  workflowSelect.value = state.draft.workflow;

  const runtimeSelect = el(
    "select",
    {
      id: "runtime",
      onchange: (e) => {
        state.draft.runtime = e.target.value;
        persistDraft();
        rerender();
      },
    },
    allowedRuntimes.map((runtime) =>
      el("option", {
        value: runtime,
        text: runtime === "apptainer" ? "Apptainer" : runtime[0].toUpperCase() + runtime.slice(1),
      })
    )
  );
  runtimeSelect.value = state.draft.runtime;

  const containerFieldsVisible = state.draft.runtime === "docker" || state.draft.runtime === "apptainer";

  form.appendChild(
    el("div", { class: "card" }, [
      el("div", { class: "card-header" }, [
        el("div", {}, [
          el("h3", { class: "card-title", text: "Submit a new run" }),
          el("div", { class: "muted", text: "Create a run, stage config overrides, and launch Snakemake." }),
        ]),
        el("div", { class: "row" }, [
          el(
            "button",
            {
              class: "btn btn-secondary",
              type: "button",
              onclick: () => applyBundledTestRunDefaults(),
            },
            "Use bundled test run"
          ),
          el("button", { class: "btn btn-primary", type: "submit" }, "Submit Run"),
        ]),
      ]),
      el("div", { class: "card-body" }, [
        el("div", { class: "help" }, [
          "Quick start: click ",
          el("span", { class: "kbd", text: "Use bundled test run" }),
          " to prefill Docker + bundled test FASTQs (",
          el("span", { class: "kbd", text: BUNDLED_TEST.inputDir }),
          ").",
        ]),
        el("div", { class: "hr" }),
        el("div", { class: "grid grid-2" }, [
          inputField({
            label: "Preset",
            help: "Prefills a minimal set of config overrides. You can refine them in the Config tab.",
            input: presetSelect,
          }),
          inputField({
            label: "Workflow",
            help:
              "If OTU is not installed on this server, submission may fail (missing Snakefile).",
            input: workflowSelect,
          }),
          inputField({
            label: "Run name",
            help: "A short label to recognize this run later.",
            input: el("input", {
              id: "run_name",
              type: "text",
              value: state.draft.run_name,
              placeholder: "e.g. COI batch 01",
              oninput: (e) => {
                state.draft.run_name = e.target.value;
                persistDraft();
              },
              required: true,
            }),
          }),
          inputField({
            label: "Runtime",
            help: "Runs are executed in containerized mode configured by the server.",
            input: runtimeSelect,
          }),
          inputField({
            label: "Input directory",
            help: "Path to a folder containing FASTQ files (must be reachable by the control node).",
            input: el("input", {
              id: "input_dir",
              type: "text",
              value: state.draft.input_dir,
              placeholder: "e.g. /workspace/tests/testing_data",
              oninput: (e) => {
                state.draft.input_dir = e.target.value;
                persistDraft();
              },
              required: true,
            }),
          }),
          inputField({
            label: "Sample source",
            help: "Folder auto-detects samples; CSV lets you define samples manually.",
            input: (() => {
              const sel = el(
                "select",
                {
                  id: "sample_source",
                  onchange: (e) => {
                    state.draft.sample_source = e.target.value;
                    persistDraft();
                    rerender();
                  },
                },
                [
                  el("option", { value: "folder", text: "Folder" }),
                  el("option", { value: "csv", text: "CSV" }),
                ]
              );
              sel.value = state.draft.sample_source;
              return sel;
            })(),
          }),
          state.draft.sample_source === "csv"
            ? inputField({
                label: "Samples CSV path",
                help: "Path to a CSV file accessible to the control node.",
                input: el("input", {
                  id: "samples_csv",
                  type: "text",
                  value: state.draft.samples_csv,
                  placeholder: "e.g. /workspace/samples.csv",
                  oninput: (e) => {
                    state.draft.samples_csv = e.target.value;
                    persistDraft();
                  },
                }),
              })
            : null,
        ]),

        el("details", { open: null }, [
          el("summary", { text: "Modules" }),
          el("div", { class: "details-body" }, [
            el("div", { class: "help" }, [
              "These toggles set ",
              el("span", { class: "kbd", text: "config_overrides.modules.*" }),
              " without affecting your other config overrides.",
            ]),
            el("div", { class: "hr" }),
            ...MODULE_DEFS.map((m) => {
              const id = `mod_${m.key}`;
              const checkbox = el("input", {
                id,
                type: "checkbox",
                checked: effectiveModuleValue(state.draft.workflow, m.key) ? true : null,
                onchange: (e) => {
                  setModuleOverride(state.draft.workflow, m.key, Boolean(e.target.checked));
                },
              });
              return el("div", { class: "field" }, [
                el("div", { class: "row" }, [checkbox, el("label", { for: id, text: m.label })]),
                el("div", { class: "help" }, m.help),
              ]);
            }),
          ]),
        ]),

        el("details", { open: containerFieldsVisible ? true : null }, [
          el("summary", { text: "Container runtime options" }),
          el("div", { class: "details-body" }, [
            el("div", { class: "grid grid-2" }, [
              inputField({
                label: "Container image",
                help: "Required for Docker/Apptainer, e.g. docker://metaworks:latest or /path/to/image.sif.",
                input: el("input", {
                  id: "container_image",
                  type: "text",
                  value: state.draft.container_image,
                  placeholder: "docker://metaworks:latest",
                  oninput: (e) => {
                    state.draft.container_image = e.target.value;
                    persistDraft();
                  },
                  disabled: !containerFieldsVisible,
                }),
              }),
              inputField({
                label: "Cache / prefix directory",
                help: "Optional cache directory for container pulls (Apptainer/Singularity).",
                input: el("input", {
                  id: "cache_dir",
                  type: "text",
                  value: state.draft.cache_dir,
                  placeholder: "e.g. /workspace/runtime/cache",
                  oninput: (e) => {
                    state.draft.cache_dir = e.target.value;
                    persistDraft();
                  },
                  disabled: !containerFieldsVisible,
                }),
              }),
            ]),
            inputField({
              label: "Bind paths (one per line)",
              help: "Optional extra bind mounts in src:dest format (Docker -v or Apptainer -B).",
              input: (() => {
                const area = el("textarea", {
                  id: "bind_paths",
                  placeholder: "/data:/data\n/refs:/refs",
                  oninput: (e) => {
                    state.draft.bind_paths_text = e.target.value;
                    persistDraft();
                  },
                  disabled: !containerFieldsVisible,
                });
                area.value = state.draft.bind_paths_text || "";
                return area;
              })(),
            }),
          ]),
        ]),

        el("details", {}, [
          el("summary", { text: "Resources and notes" }),
          el("div", { class: "details-body" }, [
            el("div", { class: "grid grid-2" }, [
              inputField({
                label: "CPU cores (optional)",
                help: "Overrides server defaults for this run.",
                input: el("input", {
                  id: "cores",
                  type: "number",
                  min: "1",
                  value: state.draft.cores,
                  placeholder: "e.g. 16",
                  oninput: (e) => {
                    state.draft.cores = e.target.value;
                    persistDraft();
                  },
                }),
              }),
              inputField({
                label: "Memory (GB, optional)",
                help: "Overrides server defaults for this run.",
                input: el("input", {
                  id: "mem_gb",
                  type: "number",
                  min: "1",
                  value: state.draft.mem_gb,
                  placeholder: "e.g. 120",
                  oninput: (e) => {
                    state.draft.mem_gb = e.target.value;
                    persistDraft();
                  },
                }),
              }),
            ]),
            inputField({
              label: "Notes (optional)",
              help: "Stored with the run metadata on the control node.",
              input: el("input", {
                id: "notes",
                type: "text",
                value: state.draft.notes,
                placeholder: "e.g. replicate run with updated classifier",
                oninput: (e) => {
                  state.draft.notes = e.target.value;
                  persistDraft();
                },
              }),
            }),
            el("div", { class: "row" }, [
              el("label", { for: "dry_run" }, "Dry-run (stage only)"),
              el("input", {
                id: "dry_run",
                type: "checkbox",
                checked: state.draft.dry_run ? true : null,
                onchange: (e) => {
                  state.draft.dry_run = Boolean(e.target.checked);
                  persistDraft();
                },
              }),
              el("div", { class: "help" }, [
                "Creates the run directory and rendered config but does not launch Snakemake.",
              ]),
            ]),
            el("div", { class: "help" }, [
              retentionPolicy === "immediate"
                ? "Server retention policy is immediate cleanup after completion."
                : retentionPolicy === "manual"
                  ? "Server retention policy is manual cleanup (files persist until deleted)."
                  : "Run files are retained until you click Download artifacts in the Runs tab, then cleaned automatically.",
            ]),
          ]),
        ]),

        el("div", { class: "hr" }),
        el("div", { class: "row" }, [
          el("a", { class: "btn btn-accent", href: "#/config" }, "Edit Config Overrides"),
          el("a", { class: "btn btn-secondary", href: "#/runs" }, "Go to Runs"),
        ]),
      ]),
    ])
  );

  root.appendChild(
    renderPageHeader(
      "Submit Run",
      "Create a pipeline run using Conda, Docker, or Apptainer. Use presets, then refine config overrides."
    )
  );
  root.appendChild(form);
  root.appendChild(responseCard);
  return root;
}

function renderRunsView() {
  const root = el("div", { class: "grid", role: "region", "aria-label": "Runs view" });

  const addInput = el("input", { type: "text", id: "add_run_id", placeholder: "esv-YYYYMMDDhhmmss" });
  const addBtn = el(
    "button",
    {
      class: "btn btn-primary",
      type: "button",
      onclick: () => {
        const id = String(addInput.value || "").trim();
        if (!id) return;
        state.trackedRuns.unshift(id);
        addInput.value = "";
        persistTrackedRuns();
        render();
      },
    },
    "Add"
  );

  const refreshBtn = el(
    "button",
    {
      class: "btn btn-secondary",
      type: "button",
      onclick: () => render(),
    },
    "Refresh"
  );

  const top = el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h3", { class: "card-title", text: "Tracked runs" }),
        el("div", { class: "muted", text: "Run IDs are persisted in your browser (localStorage)." }),
      ]),
      el("div", { class: "row" }, [refreshBtn]),
    ]),
    el("div", { class: "card-body" }, [
      el("div", { class: "row" }, [
        el("div", { class: "field", style: "flex: 1; min-width: 240px; margin: 0" }, [
          el("label", { for: "add_run_id", text: "Add run by ID" }),
          addInput,
        ]),
        addBtn,
        el(
          "button",
          {
            class: "btn btn-danger",
            type: "button",
            onclick: () => {
              if (!confirm("Clear all tracked runs from this browser?")) return;
              state.trackedRuns = [];
              persistTrackedRuns();
              render();
            },
          },
          "Clear"
        ),
      ]),
    ]),
  ]);

  const list = el("div", { class: "grid" });

  const runIds = state.trackedRuns;
  if (runIds.length === 0) {
    list.appendChild(
      el("div", { class: "card" }, [
        el("div", { class: "card-body" }, [
          el("div", { class: "muted" }, [
            "No tracked runs yet. Submit a run, or paste a ",
            el("span", { class: "kbd", text: "run_id" }),
            " here to track it.",
          ]),
        ]),
      ])
    );
  }

  for (const runId of runIds) {
    const card = el("div", { class: "card" });
    const headerRight = el("div", { class: "row" }, [el("span", { class: "pill neutral", text: "loading" })]);
    const header = el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h3", { class: "card-title", text: runId }),
        el("div", { class: "muted", text: "Fetching status…" }),
      ]),
      headerRight,
    ]);
    const body = el("div", { class: "card-body" }, [el("div", { class: "muted", text: "Loading…" })]);
    card.appendChild(header);
    card.appendChild(body);
    list.appendChild(card);

    (async () => {
      try {
        const status = await api.getRunStatus(runId);
        state.runStatusById = state.runStatusById || {};
        state.runStatusById[runId] = status;
        const prog = progressForStatus(status.status);

        header.replaceWith(
          el("div", { class: "card-header" }, [
            el("div", {}, [
              el("h3", { class: "card-title", text: status.run_name ? `${status.run_name} · ${runId}` : runId }),
              el(
                "div",
                { class: "muted" },
                `Workflow: ${status.workflow || "—"} · Runtime: ${status.runtime || "—"}`
              ),
            ]),
            el("div", { class: "row" }, [statusPill(status.status)]),
          ])
        );

        const progress = el("div", { class: `progress ${prog.kind === "indeterminate" ? "indeterminate" : ""}` }, [
          el("div", { style: prog.kind === "determinate" ? `width:${prog.pct}%` : "" }),
        ]);

        const logsPre = el("pre", { class: "code", text: "" });
        const logsBtn = el(
          "button",
          {
            class: "btn btn-secondary",
            type: "button",
            onclick: async () => {
              logsBtn.disabled = true;
              logsBtn.textContent = "Loading logs…";
              try {
                const logs = await api.getRunLogs(runId);
                logsPre.textContent = Array.isArray(logs.tail) ? logs.tail.join("\n") : String(logs.tail || "");
              } catch (err) {
                logsPre.textContent = `Error: ${String(err.message || err)}`;
              } finally {
                logsBtn.disabled = false;
                logsBtn.textContent = "Tail logs";
              }
            },
          },
          "Tail logs"
        );

        const shouldAutoLoadLogs = ["failed", "cancelled", "completed"].includes(
          String(status.status || "").toLowerCase()
        );
        if (shouldAutoLoadLogs) {
          (async () => {
            try {
              const logs = await api.getRunLogs(runId);
              logsPre.textContent = Array.isArray(logs.tail)
                ? logs.tail.join("\n")
                : String(logs.tail || "");
            } catch (err) {
              logsPre.textContent = `Error: ${String(err.message || err)}`;
            }
          })();
        }

        const cancelBtn = el(
          "button",
          {
            class: "btn btn-danger",
            type: "button",
            disabled: String(status.status || "").toLowerCase() !== "running" ? true : null,
            onclick: async () => {
              if (!confirm(`Cancel run ${runId}?`)) return;
              try {
                await api.cancelRun(runId);
                toast("Run cancelled", runId);
                render();
              } catch (err) {
                toast("Cancel failed", String(err.message || err));
              }
            },
          },
          "Cancel"
        );

        const removeBtn = el(
          "button",
          {
            class: "btn btn-secondary",
            type: "button",
            onclick: () => {
              state.trackedRuns = state.trackedRuns.filter((id) => id !== runId);
              persistTrackedRuns();
              render();
            },
          },
          "Remove"
        );

        const deleteBtn = el(
          "button",
          {
            class: "btn btn-danger",
            type: "button",
            onclick: async () => {
              if (!confirm(`Delete run record and files for ${runId}?`)) return;
              try {
                await api.deleteRun(runId);
                state.trackedRuns = state.trackedRuns.filter((id) => id !== runId);
                persistTrackedRuns();
                toast("Run deleted", runId);
                render();
              } catch (err) {
                toast("Delete failed", String(err.message || err));
              }
            },
          },
          "Delete"
        );

        const downloadLogBtn = el(
          "a",
          {
            class: "btn btn-accent",
            href: `/runs/${encodeURIComponent(runId)}/logs/download`,
          },
          "Download log"
        );
        const downloadArtifactsBtn = el(
          "a",
          {
            class: "btn btn-accent",
            href: `/runs/${encodeURIComponent(runId)}/artifacts/download`,
          },
          "Download artifacts"
        );
        const actions = el("div", { class: "row" }, [
          logsBtn,
          ...(status.keep_outputs ? [downloadLogBtn, downloadArtifactsBtn] : []),
          cancelBtn,
          removeBtn,
          deleteBtn,
        ]);

        const metaGrid = el("div", { class: "grid grid-2" }, [
          el("div", { class: "field" }, [
            el("label", { text: "Submitted" }),
            el("div", { class: "help", text: formatTime(status.submitted_at) }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Started" }),
            el("div", { class: "help", text: formatTime(status.started_at) }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Ended" }),
            el("div", { class: "help", text: formatTime(status.ended_at) }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Message" }),
            el("div", { class: "help", text: status.message || "—" }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Return code" }),
            el("div", { class: "help", text: status.return_code == null ? "—" : String(status.return_code) }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "PID" }),
            el("div", { class: "help", text: status.pid == null ? "—" : String(status.pid) }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Run dir" }),
            el("div", { class: "help", text: status.run_dir || "—" }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Log path" }),
            el("div", { class: "help", text: status.log_path || "—" }),
          ]),
          el("div", { class: "field" }, [
            el("label", { text: "Persistence" }),
            el(
              "div",
              {
                class: "help",
                text: status.keep_outputs ? "Retained (download/debug enabled)" : "Ephemeral (auto-cleanup)",
              }
            ),
          ]),
        ]);

        const cmdDetails =
          status.command
            ? el("details", {}, [
                el("summary", { text: "Command" }),
                el("div", { class: "details-body" }, [el("pre", { class: "code", text: status.command })]),
              ])
            : null;

        const logsDetails = el(
          "details",
          {
            open:
              String(status.status || "").toLowerCase() === "failed"
                ? true
                : null,
          },
          [
          el("summary", { text: "Log tail" }),
          el("div", { class: "details-body" }, [logsPre]),
          ]
        );

        body.replaceWith(
          el("div", { class: "card-body" }, [
            progress,
            el("div", { class: "hr" }),
            metaGrid,
            cmdDetails,
            logsDetails,
            el("div", { class: "hr" }),
            actions,
          ])
        );
      } catch (err) {
        state.runStatusById = state.runStatusById || {};
        delete state.runStatusById[runId];
        body.replaceWith(
          el("div", { class: "card-body" }, [
            el("div", { class: "row" }, [statusPill("unknown")]),
            el("div", { class: "help", text: `Failed to fetch status: ${String(err.message || err)}` }),
            el("div", { class: "hr" }),
            el("div", { class: "row" }, [
              el(
                "button",
                {
                  class: "btn btn-secondary",
                  type: "button",
                  onclick: () => render(),
                },
                "Retry"
              ),
              el(
                "button",
                {
                  class: "btn btn-secondary",
                  type: "button",
                  onclick: () => {
                    state.trackedRuns = state.trackedRuns.filter((id) => id !== runId);
                    persistTrackedRuns();
                    render();
                  },
                },
                "Remove"
              ),
            ]),
          ])
        );
      }
    })();
  }

  root.appendChild(
    renderPageHeader(
      "Runs",
      "Track active/completed runs by ID. Active runs are polled every 5s without full page refresh."
    )
  );
  root.appendChild(top);
  root.appendChild(list);

  // Poll active runs while viewing Runs without rebuilding the full page each cycle.
  if (state.view.intervalId) clearInterval(state.view.intervalId);
  state.view.intervalId = setInterval(() => {
    pollActiveRunsAndMaybeRender().catch(() => {});
  }, 5000);

  return root;
}

function renderConfigView() {
  const root = el("div", { class: "grid", role: "region", "aria-label": "Config view" });

  const workflow = state.draft.workflow;
  const overrides = currentOverrides();

  const headerRight = el("div", { class: "row" }, [
    el(
      "button",
      {
        class: "btn btn-secondary",
        type: "button",
        onclick: () => {
          state.overridesByWorkflow[workflow] = applyPreset(state.draft.preset, workflow);
          persistOverrides(workflow);
          toast("Config", "Preset overrides applied.");
          render();
        },
      },
      "Reapply preset"
    ),
    el(
      "button",
      {
        class: "btn btn-danger",
        type: "button",
        onclick: () => {
          if (!confirm(`Clear overrides for ${workflow.toUpperCase()}?`)) return;
          state.overridesByWorkflow[workflow] = {};
          persistOverrides(workflow);
          toast("Config", "Overrides cleared.");
          render();
        },
      },
      "Reset overrides"
    ),
  ]);

  const defaultsCard = el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h3", { class: "card-title", text: `Defaults (${workflow.toUpperCase()})` }),
        el("div", { class: "muted", text: "Loaded from the backend and used to compute minimal overrides." }),
      ]),
      el("div", { class: "row" }, [
        el(
          "button",
          {
            class: "btn btn-primary",
            type: "button",
            onclick: async () => {
              try {
                await ensureDefaults(workflow);
                toast("Config", "Defaults loaded.");
                render();
              } catch (err) {
                toast("Config", `Failed to load defaults: ${String(err.message || err)}`);
              }
            },
          },
          "Load defaults"
        ),
      ]),
    ]),
    el("div", { class: "card-body" }, [
      el("div", { class: "help" }, [
        "Tip: Use ",
        el("span", { class: "kbd", text: "Load defaults" }),
        " to edit values with a structured form, or edit overrides JSON directly below.",
      ]),
    ]),
  ]);

  const overridesText = el("textarea", { id: "overrides_json", spellcheck: "false" });
  overridesText.value = JSON.stringify(overrides, null, 2);

  const overridesCard = el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h3", { class: "card-title", text: "Overrides (JSON)" }),
        el("div", { class: "muted", text: "Only changed values are sent to the backend." }),
      ]),
      el("div", { class: "row" }, [
        el(
          "button",
          {
            class: "btn btn-secondary",
            type: "button",
            onclick: async () => {
              try {
                await navigator.clipboard.writeText(overridesText.value);
                toast("Clipboard", "Overrides copied.");
              } catch (err) {
                toast("Clipboard", `Copy failed: ${String(err.message || err)}`);
              }
            },
          },
          "Copy"
        ),
        el(
          "button",
          {
            class: "btn btn-primary",
            type: "button",
            onclick: () => {
              const parsed = safeJsonParse(overridesText.value);
              if (!parsed.ok) {
                toast("Invalid JSON", String(parsed.error.message || parsed.error));
                return;
              }
              state.overridesByWorkflow[workflow] = parsed.value || {};
              persistOverrides(workflow);
              toast("Config", "Overrides updated.");
              render();
            },
          },
          "Apply"
        ),
      ]),
    ]),
    el("div", { class: "card-body" }, [overridesText]),
  ]);

  const renderedOut = el("pre", { class: "code", text: "" });
  const renderBtn = el(
    "button",
    {
      class: "btn btn-accent",
      type: "button",
      onclick: async () => {
        renderBtn.disabled = true;
        renderBtn.textContent = "Rendering…";
        try {
          const res = await api.renderConfig(workflow, currentOverrides());
          renderedOut.textContent = res.merged || "";
        } catch (err) {
          renderedOut.textContent = `Error: ${String(err.message || err)}`;
        } finally {
          renderBtn.disabled = false;
          renderBtn.textContent = "Preview rendered config";
        }
      },
    },
    "Preview rendered config"
  );

  const previewCard = el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h3", { class: "card-title", text: "Rendered config preview" }),
        el("div", { class: "muted", text: "Backend merge result (defaults + overrides)." }),
      ]),
      el("div", { class: "row" }, [renderBtn]),
    ]),
    el("div", { class: "card-body" }, [renderedOut]),
  ]);

  const defaults = state.defaultsByWorkflow[workflow];
  const formCard = el("div", { class: "card" }, [
    el("div", { class: "card-header" }, [
      el("div", {}, [
        el("h3", { class: "card-title", text: "Structured editor" }),
        el("div", { class: "muted", text: "Edits values against the loaded defaults." }),
      ]),
      el("div", { class: "row" }, [
        el("a", { class: "btn btn-secondary", href: "#/submit" }, "Back to Submit"),
      ]),
    ]),
    el("div", { class: "card-body" }, [
      defaults
        ? renderConfigForm(defaults, workflow)
        : el("div", { class: "help" }, [
            "Load defaults to use the structured editor. If defaults are missing on the server, you can still edit overrides JSON above.",
          ]),
    ]),
  ]);

  root.appendChild(renderPageHeader("Config", "Edit workflow defaults and produce minimal JSON overrides.", headerRight));
  root.appendChild(defaultsCard);
  root.appendChild(formCard);
  root.appendChild(overridesCard);
  root.appendChild(previewCard);

  return root;
}

function renderConfigForm(defaults, workflow) {
  const overrides = state.overridesByWorkflow[workflow] || {};
  const wrapper = el("div", { class: "grid" });

  function addField(path, defaultValue) {
    const key = path[path.length - 1];
    const overrideValue = getAt(overrides, path);
    const value = mergeViewValue(defaultValue, overrideValue);

    const id = `cfg_${path.join("__")}`;

    if (typeof defaultValue === "boolean") {
      const input = el("input", {
        id,
        type: "checkbox",
        checked: value ? true : null,
        onchange: (e) => {
          const next = Boolean(e.target.checked);
          if (deepEqual(next, defaultValue)) deleteAt(overrides, path);
          else setAt(overrides, path, next);
          state.overridesByWorkflow[workflow] = overrides;
          persistOverrides(workflow);
          syncOverridesTextarea();
        },
      });
      return el("div", { class: "field" }, [
        el("label", { for: id, text: key }),
        el("div", { class: "row" }, [input, el("div", { class: "help", text: path.join(".") })]),
      ]);
    }

    if (typeof defaultValue === "number") {
      const input = el("input", {
        id,
        type: "number",
        value: String(value),
        oninput: (e) => {
          const raw = e.target.value;
          const next = raw === "" ? defaultValue : Number(raw);
          if (Number.isNaN(next)) return;
          if (deepEqual(next, defaultValue)) deleteAt(overrides, path);
          else setAt(overrides, path, next);
          state.overridesByWorkflow[workflow] = overrides;
          persistOverrides(workflow);
          syncOverridesTextarea();
        },
      });
      return el("div", { class: "field" }, [
        el("label", { for: id, text: key }),
        input,
        el("div", { class: "help", text: path.join(".") }),
      ]);
    }

    if (Array.isArray(defaultValue)) {
      const input = el("textarea", {
        id,
        oninput: (e) => {
          const lines = String(e.target.value || "")
            .split("\n")
            .map((s) => s.trim())
            .filter((s) => s.length > 0);
          const next = lines;
          if (deepEqual(next, defaultValue)) deleteAt(overrides, path);
          else setAt(overrides, path, next);
          state.overridesByWorkflow[workflow] = overrides;
          persistOverrides(workflow);
          syncOverridesTextarea();
        },
      });
      input.value = Array.isArray(value) ? value.join("\n") : "";
      return el("div", { class: "field" }, [
        el("label", { for: id, text: key }),
        input,
        el("div", { class: "help", text: `${path.join(".")} (one item per line)` }),
      ]);
    }

    const input = el("input", {
      id,
      type: "text",
      value: value == null ? "" : String(value),
      oninput: (e) => {
        const next = e.target.value;
        if (deepEqual(next, defaultValue == null ? "" : String(defaultValue))) deleteAt(overrides, path);
        else setAt(overrides, path, next);
        state.overridesByWorkflow[workflow] = overrides;
        persistOverrides(workflow);
        syncOverridesTextarea();
      },
    });
    return el("div", { class: "field" }, [
      el("label", { for: id, text: key }),
      input,
      el("div", { class: "help", text: path.join(".") }),
    ]);
  }

  function addGroup(path, obj) {
    const title = path.length === 0 ? "root" : path[path.length - 1];
    const details = el("details", { open: path.length < 2 ? true : null }, [
      el("summary", { text: title }),
      el("div", { class: "details-body" }, []),
    ]);
    const body = details.querySelector(".details-body");
    const keys = Object.keys(obj || {});
    keys.sort();
    for (const k of keys) {
      const child = obj[k];
      const nextPath = path.concat([k]);
      if (child && typeof child === "object" && !Array.isArray(child)) {
        body.appendChild(addGroup(nextPath, child));
      } else {
        body.appendChild(addField(nextPath, child));
      }
    }
    return details;
  }

  wrapper.appendChild(addGroup([], defaults));
  return wrapper;
}

function renderAssetsView() {
  const root = el("div", { class: "grid", role: "region", "aria-label": "Assets view" });

  function makeAssetSection({ title, listFn, uploadPath, deleteFn, hint }) {
    const listArea = el("div", { class: "grid" });
    const fileInput = el("input", { type: "file", multiple: true, style: "display:none" });
    const drop = el("div", { class: "dropzone", tabindex: "0" }, [
      el("div", { class: "row" }, [
        el("span", { class: "pill neutral", text: title }),
        el("span", { class: "muted", text: hint }),
      ]),
      el("div", { class: "row" }, [
        el(
          "button",
          {
            class: "btn btn-primary",
            type: "button",
            onclick: () => fileInput.click(),
          },
          "Browse files"
        ),
        el(
          "button",
          {
            class: "btn btn-secondary",
            type: "button",
            onclick: () => refresh(),
          },
          "Refresh"
        ),
      ]),
      fileInput,
    ]);

    async function refresh() {
      listArea.innerHTML = "";
      try {
        const res = await listFn();
        const items = (res && res.items) || [];
        if (!items.length) {
          listArea.appendChild(el("div", { class: "muted", text: "No files uploaded." }));
          return;
        }
        for (const name of items.sort()) {
          listArea.appendChild(
            el("div", { class: "card" }, [
              el("div", { class: "card-body" }, [
                el("div", { class: "row" }, [
                  el("div", { style: "flex:1; min-width: 220px" }, [
                    el("div", { style: "font-weight: 800" }, name),
                    el("div", { class: "muted", text: "Stored in runtime directory on the control node." }),
                  ]),
                  el(
                    "button",
                    {
                      class: "btn btn-danger",
                      type: "button",
                      onclick: async () => {
                        if (!confirm(`Delete ${name}?`)) return;
                        try {
                          await deleteFn(name);
                          toast("Deleted", name);
                          refresh();
                        } catch (err) {
                          toast("Delete failed", String(err.message || err));
                        }
                      },
                    },
                    "Delete"
                  ),
                ]),
              ]),
            ])
          );
        }
      } catch (err) {
        listArea.appendChild(el("div", { class: "muted", text: `Failed to list: ${String(err.message || err)}` }));
      }
    }

    function uploadFiles(files) {
      for (const file of Array.from(files || [])) {
        const row = el("div", { class: "card" }, [
          el("div", { class: "card-body" }, [
            el("div", { style: "font-weight: 800" }, file.name),
            el("div", { class: "progress" }, [el("div", { style: "width:0%" })]),
            el("div", { class: "muted", text: "Uploading…" }),
          ]),
        ]);
        listArea.prepend(row);
        const bar = row.querySelector(".progress > div");
        const status = row.querySelector(".muted");

        const xhr = new XMLHttpRequest();
        xhr.open("POST", `${API_BASE}${uploadPath}`);
        xhr.upload.onprogress = (evt) => {
          if (!evt.lengthComputable) return;
          const pct = Math.round((evt.loaded / evt.total) * 100);
          bar.style.width = `${pct}%`;
          status.textContent = `Uploading… ${pct}%`;
        };
        xhr.onload = () => {
          if (xhr.status >= 200 && xhr.status < 300) {
            bar.style.width = "100%";
            status.textContent = "Uploaded";
            toast("Uploaded", file.name);
            refresh();
          } else {
            status.textContent = `Error: HTTP ${xhr.status}`;
            toast("Upload failed", `${file.name}: HTTP ${xhr.status}`);
          }
        };
        xhr.onerror = () => {
          status.textContent = "Error";
          toast("Upload failed", `${file.name}: network error`);
        };
        const data = new FormData();
        data.append("file", file);
        xhr.send(data);
      }
    }

    fileInput.addEventListener("change", (e) => uploadFiles(e.target.files));
    drop.addEventListener("dragover", (e) => {
      e.preventDefault();
      drop.classList.add("dragover");
    });
    drop.addEventListener("dragleave", () => drop.classList.remove("dragover"));
    drop.addEventListener("drop", (e) => {
      e.preventDefault();
      drop.classList.remove("dragover");
      uploadFiles(e.dataTransfer.files);
    });
    drop.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") fileInput.click();
    });

    const section = el("div", { class: "card" }, [
      el("div", { class: "card-header" }, [
        el("div", {}, [el("h3", { class: "card-title", text: title }), el("div", { class: "muted", text: hint })]),
        el("div"),
      ]),
      el("div", { class: "card-body" }, [drop, el("div", { class: "hr" }), listArea]),
    ]);

    refresh();
    return section;
  }

  root.appendChild(renderPageHeader("Assets", "Upload and manage adapters and classifiers used by the pipeline."));
  root.appendChild(
    makeAssetSection({
      title: "Adapters",
      listFn: api.listAdapters,
      uploadPath: "/adapters",
      deleteFn: api.deleteAdapter,
      hint: "Upload adapter FASTA files (Cutadapt).",
    })
  );
  root.appendChild(
    makeAssetSection({
      title: "Classifiers",
      listFn: api.listClassifiers,
      uploadPath: "/classifiers",
      deleteFn: api.deleteClassifier,
      hint: "Upload classifier property files (e.g. RDP).",
    })
  );
  return root;
}

function routeFromHash() {
  const hash = window.location.hash || "#/submit";
  const parts = hash.replace(/^#\//, "").split("/");
  const route = parts[0] || "submit";
  if (!["submit", "runs", "config", "assets"].includes(route)) return "submit";
  return route;
}

function updateNav(route) {
  for (const link of document.querySelectorAll(".nav-link")) {
    const r = link.getAttribute("data-route");
    if (r === route) link.setAttribute("aria-current", "page");
    else link.removeAttribute("aria-current");
  }
}

function render() {
  const route = routeFromHash();
  state.view.route = route;
  updateNav(route);

  const appRoot = document.getElementById("app");
  appRoot.innerHTML = "";

  if (route === "submit") appRoot.appendChild(renderSubmitView());
  else if (route === "runs") appRoot.appendChild(renderRunsView());
  else if (route === "config") appRoot.appendChild(renderConfigView());
  else if (route === "assets") appRoot.appendChild(renderAssetsView());
}

window.addEventListener("hashchange", () => {
  if (state.view.intervalId) {
    clearInterval(state.view.intervalId);
    state.view.intervalId = null;
  }
  render();
});

// First render
persistTrackedRuns();
if (!state.overridesByWorkflow.esv || Object.keys(state.overridesByWorkflow.esv).length === 0) {
  state.overridesByWorkflow.esv = applyPreset(state.draft.preset || "coi", "esv");
  persistOverrides("esv");
}
if (!state.overridesByWorkflow.otu || Object.keys(state.overridesByWorkflow.otu).length === 0) {
  state.overridesByWorkflow.otu = applyPreset(state.draft.preset || "coi", "otu");
  persistOverrides("otu");
}
normalizeDraftRuntime();

(async () => {
  try {
    const settingsRes = await api.getSettingsPaths();
    const allowed = Array.isArray(settingsRes?.allowed_runtimes)
      ? settingsRes.allowed_runtimes.map((v) => String(v).toLowerCase()).filter(Boolean)
      : [];
    if (allowed.length) state.systemSettings.allowed_runtimes = allowed;
    state.systemSettings.retention_policy = String(
      settingsRes?.retention_policy || state.systemSettings.retention_policy
    ).toLowerCase();
    normalizeDraftRuntime();
  } catch {
    // Keep local defaults if settings endpoint is unavailable.
  } finally {
    render();
  }
})();
