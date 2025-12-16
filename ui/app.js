const apiBase = window.location.origin;
document.getElementById("api-base").innerText = apiBase;

const runtimeSelect = document.getElementById("runtime-select");
const containerFields = document.getElementById("container-fields");
const cacheDirInput = document.getElementById("cache-dir");
const inputDirInput = document.getElementById("input-dir");
const attachRunInput = document.getElementById("attach-run-input");
const attachRunBtn = document.getElementById("attach-run-btn");
const runsContainer = document.getElementById("runs-container");
const refreshAllBtn = document.getElementById("refresh-all-btn");
const pollingCount = document.getElementById("polling-count");
const emptyRuns = document.getElementById("empty-runs");
const tabs = document.querySelectorAll(".tab");
const views = document.querySelectorAll(".view");
const runForm = document.getElementById("run-form");
const workflowSelect = runForm.querySelector("select[name='workflow']");
const configCardsEl = document.getElementById("config-cards");
const CONFIG_LAYOUT = [
  { id: "SEQPREP", label: "SeqPrep", fields: ["SEQPREP"] },
  { id: "CUTADAPT", label: "Primer Trimming", fields: ["CUTADAPT"] },
  { id: "DENOISING", label: "Denoising", fields: ["pooling", "VSEARCH_DENOISE", "VSEARCH_TABLE"] },
  { id: "ITSX", label: "ITSx extractor", fields: ["ITSpart"], hint: "Skip if not needed; leave defaults." },
  { id: "TAXONOMY", label: "Taxonomic assignment", fields: ["marker", "RDP"] },
  {
    id: "PSEUDOGENE",
    label: "Pseudogene filtering",
    fields: ["pseudogene_filtering", "grep_type", "taxon1", "taxon2", "removal_type", "hmm", "ORFFINDER"],
  },
  { id: "REPORT", label: "Report type", fields: ["report_type"] },
];
const CONFIG_DOCS = {
  SEQPREP: {
    q: { hint: "Phred quality cutoff." },
    o: { hint: "Minimum overlap length between reads." },
    m: { hint: "Maximum mismatch fraction in overlap." },
    n: { hint: "Minimum matching overlap fraction." },
  },
  CUTADAPT: {
    fasta: { hint: "FASTA file with linked adapters." },
    m: { hint: "Minimum sequence length to retain after trimming." },
    q: { hint: "Phred quality cutoffs at ends (e.g., 20,20)." },
    e: { hint: "Error rate." },
    O: { hint: "Minimum adapter overlap." },
    mn: { hint: "Max number of Ns." },
    rc: { hint: "Enable reverse complement search.", options: ["Yes", "No"] },
  },
  pooling: { hint: "Pool samples for denoising.", options: ["Yes", "No"] },
  VSEARCH_DENOISE: {
    minsize: { hint: "Minimum reads per cluster to retain." },
  },
  VSEARCH_TABLE: {
    t: { hint: "Threads to use when building the table." },
  },
  marker: { hint: "Marker classifier to use." },
  ITSpart: { hint: "ITS spacer region; leave as-is if not using ITS." },
  RDP: {
    memory: { hint: "Memory to allocate to RDP classifier (e.g., 20g)." },
    custom: { hint: "Use custom-trained dataset?", options: ["yes", "no"] },
    t: { hint: "Path to custom rRNAClassifier.properties file." },
    c: { hint: "Confidence threshold (16S built-in)." },
    f: { hint: "Output format flag (16S built-in)." },
    g: { hint: "Built-in fungal classifier choice." },
  },
  pseudogene_filtering: { hint: "Enable pseudogene filtering?", options: ["yes", "no"] },
  grep_type: { hint: "1 = simple grep, 2 = compound grep." },
  taxon1: { hint: "First grep filter (e.g., -e Arthropoda)." },
  taxon2: { hint: "Second grep filter for compound searches." },
  removal_type: { hint: "1 = ORF length filter, 2 = HMM score filter." },
  hmm: { hint: "HMM profile path when using removal_type=2." },
  ORFFINDER: {
    g: { hint: "Genetic code (NCBI number)." },
    s: { hint: "Start codon mode." },
    ml: { hint: "Minimum ORF length." },
    n: { hint: "Ignore nested ORFs.", options: ["true", "false"] },
    strand: { hint: "Strand to translate.", options: ["both", "plus", "minus"] },
  },
  report_type: {
    hint: "1 = single combined report CSV, 2 = separate ESV outputs (auto-zip on finish).",
    inlineNote: "(1 = single results.csv, 2 = split ESV outputs, auto-zipped)",
    options: ["1", "2"],
  },
};

const REFRESH_MS = 5000;
let autoBindValue = "";
const runCards = new Map();
let sectionOverrides = {};
let defaultSections = {};

function cloneDeep(obj) {
  if (typeof structuredClone === "function") return structuredClone(obj);
  return JSON.parse(JSON.stringify(obj));
}

function deepEqual(a, b) {
  try {
    return JSON.stringify(a) === JSON.stringify(b);
  } catch (_err) {
    return a === b;
  }
}

function getDoc(section, field) {
  const secDoc = CONFIG_DOCS[section] || {};
  if (field === null) return secDoc || {};
  if (typeof secDoc === "object" && !Array.isArray(secDoc)) {
    // Prefer nested field docs, otherwise fall back to section-level docs (for flat fields like report_type).
    if (secDoc[field]) return secDoc[field];
    if (secDoc.hint || secDoc.options) return secDoc;
  }
  return {};
}

function switchView(targetId) {
  tabs.forEach((btn) => btn.classList.toggle("active", btn.dataset.viewTarget === targetId));
  views.forEach((view) => view.classList.toggle("active", view.id === targetId));
}

tabs.forEach((btn) => {
  btn.addEventListener("click", () => switchView(btn.dataset.viewTarget));
});

function toggleContainerFields() {
  const runtime = runtimeSelect.value;
  containerFields.style.display = runtime === "conda" ? "none" : "block";
}
runtimeSelect.addEventListener("change", toggleContainerFields);
toggleContainerFields();
workflowSelect.addEventListener("change", () => loadDefaultConfig(workflowSelect.value));

runForm.addEventListener("submit", async (e) => {
  e.preventDefault();
  const data = new FormData(runForm);
  const runtime = data.get("runtime");
  const bindsText = data.get("bind_paths");
  const parsedBinds = runtime === "conda" ? [] : parseBindPaths(bindsText);
  const payload = {
    workflow: data.get("workflow"),
    run_name: data.get("run_name"),
    runtime,
    container_image: runtime === "conda" ? null : data.get("container_image") || null,
    bind_paths: runtime === "conda" ? [] : parsedBinds.length ? parsedBinds : autoBindValue ? [autoBindValue] : [],
    cache_dir: runtime === "conda" ? null : data.get("cache_dir") || null,
    input_dir: data.get("input_dir"),
    samples_csv: data.get("samples_csv") || null,
    sample_source: data.get("sample_source"),
    config_overrides: sectionOverrides,
    cores: data.get("cores") ? Number(data.get("cores")) : null,
    mem_gb: data.get("mem_gb") ? Number(data.get("mem_gb")) : null,
    dry_run: data.get("dry_run") === "on",
  };
  try {
    const res = await fetchJSON("/runs", { method: "POST", body: JSON.stringify(payload) });
    renderPre("run-response", res);
    const mergedConfig = buildCurrentConfig();
    const reportType = mergedConfig.report_type || mergedConfig.reportType || null;
    addRunToBoard(res, { reportType });
  } catch (err) {
    alert(err.message);
  }
});

attachRunBtn.addEventListener("click", async () => {
  const runId = attachRunInput.value.trim();
  if (!runId) {
    alert("Enter a run ID first");
    return;
  }
  await addRunById(runId);
});

refreshAllBtn.addEventListener("click", () => {
  runCards.forEach((_, runId) => {
    refreshRun(runId);
  });
});

function parseBindPaths(text) {
  if (!text) return [];
  return text
    .split(/[\n,]/)
    .map((v) => v.trim())
    .filter(Boolean);
}

async function loadDefaultConfig(workflow) {
  if (workflowSelect) workflowSelect.value = workflow;
  const res = await fetchJSON(`/configs/defaults/${workflow}/sections`);
  defaultSections = res.sections || {};
  sectionOverrides = {};
  renderConfigCards(defaultSections);
}

async function renderMergedConfig() {
  const data = new FormData(runForm);
  const payload = {
    workflow: data.get("workflow"),
    overrides: sectionOverrides,
  };
  const res = await fetchJSON("/configs/render", { method: "POST", body: JSON.stringify(payload) });
  renderPre("config-view", res.merged);
}

function renderConfigCards(sections) {
  configCardsEl.innerHTML = "";
  if (!Object.keys(sections || {}).length) {
    configCardsEl.innerHTML = `<p class="muted">No config sections loaded.</p>`;
    return;
  }

  CONFIG_LAYOUT.forEach((section) => {
    const presentFields = section.fields.filter((field) => sections[field] !== undefined);
    if (!presentFields.length) return;
    const card = document.createElement("div");
    card.className = "config-card";
    const header = document.createElement("div");
    header.className = "field-label";
    const hint = section.hint ? ` title="${section.hint}"` : "";
    header.innerHTML = `<span${hint}>${section.label}</span>`;
    card.appendChild(header);

    presentFields.forEach((fieldKey) => {
      const value = sections[fieldKey];
      if (value && typeof value === "object" && !Array.isArray(value)) {
        Object.entries(value).forEach(([nestedKey, nestedVal]) => {
          renderFieldRow(card, fieldKey, [nestedKey], nestedVal);
        });
      } else {
        renderFieldRow(card, fieldKey, [], value);
      }
    });
    configCardsEl.appendChild(card);
  });
}

function renderFieldRow(card, sectionKey, path, baseValue) {
  const row = document.createElement("div");
  row.className = "field-row";
  const fieldKey = path[path.length - 1] || sectionKey;
  const doc = getDoc(sectionKey, fieldKey);
  const label = document.createElement("label");
  const helper = doc.inlineNote || doc.hint || "";
  if (helper) {
    label.title = helper;
  }
  const helperTrimmed = String(helper).trim();
  const helperUnwrapped =
    helperTrimmed.startsWith("(") && helperTrimmed.endsWith(")")
      ? helperTrimmed.slice(1, -1).trim()
      : helperTrimmed;
  const hideFieldName = sectionKey === "report_type" && path.length === 0;
  label.textContent = hideFieldName && helperUnwrapped ? helperUnwrapped : fieldKey;
  if (helper && !hideFieldName) {
    const note = document.createElement("span");
    note.className = "field-note";
    const alreadyWrapped = helperTrimmed.startsWith("(") && helperTrimmed.endsWith(")");
    note.textContent = alreadyWrapped ? ` ${helperTrimmed}` : ` (${helperTrimmed})`;
    label.appendChild(note);
  }
  row.appendChild(label);

  const currentValue = getValue(sectionKey, path, baseValue);
  const options = doc.options && doc.options.length === 2 ? doc.options : null;

  if (options) {
    row.classList.add("toggle-row");
    const toggle = document.createElement("div");
    toggle.className = "toggle-group";
    options.forEach((opt) => {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.textContent = opt;
      btn.className = String(opt) === String(currentValue) ? "active" : "";
      btn.addEventListener("click", () => {
        options.forEach((_o, idx) => {
          toggle.children[idx].classList.remove("active");
        });
        btn.classList.add("active");
        updateSectionOverride(sectionKey, path, opt, baseValue);
      });
      toggle.appendChild(btn);
    });
    row.appendChild(toggle);
  } else if (typeof baseValue === "number") {
    const input = document.createElement("input");
    input.type = "number";
    input.value = currentValue ?? "";
    input.addEventListener("input", () => updateSectionOverride(sectionKey, path, input.value, baseValue));
    row.appendChild(input);
  } else if (typeof baseValue === "object") {
    const textarea = document.createElement("textarea");
    textarea.value = JSON.stringify(currentValue, null, 2);
    textarea.addEventListener("input", () => updateSectionOverride(sectionKey, path, textarea.value, baseValue, true));
    row.appendChild(textarea);
  } else {
    const input = document.createElement("input");
    input.type = "text";
    input.value = currentValue ?? "";
    input.addEventListener("input", () => updateSectionOverride(sectionKey, path, input.value, baseValue));
    row.appendChild(input);
  }
  card.appendChild(row);
}

function getValue(sourceKey, path, baseValue) {
  const overrideSection = sectionOverrides[sourceKey];
  if (!overrideSection) return baseValue;
  if (!path.length) return overrideSection;
  let cursor = overrideSection;
  for (let i = 0; i < path.length; i += 1) {
    if (cursor == null) return baseValue;
    cursor = cursor[path[i]];
  }
  return cursor === undefined ? baseValue : cursor;
}

function updateSectionOverride(sectionKey, path, rawValue, baseValue, asJson = false) {
  let parsed = rawValue;
  try {
    if (asJson) {
      parsed = JSON.parse(rawValue || "null");
    } else if (typeof baseValue === "number") {
      parsed = rawValue === "" ? null : Number(rawValue);
      if (Number.isNaN(parsed)) parsed = null;
    }
  } catch (_err) {
    return;
  }

  const baseSection = defaultSections[sectionKey];
  const isRootField = path.length === 0;

  // For scalar sections like report_type, the path is empty and we need to set the override directly.
  if (isRootField) {
    const matchesDefault = parsed === null || deepEqual(parsed, baseSection);
    if (matchesDefault) {
      delete sectionOverrides[sectionKey];
    } else {
      sectionOverrides[sectionKey] = parsed;
    }
    renderConfigCards(defaultSections);
    return;
  }

  const nextOverrides = cloneDeep(sectionOverrides[sectionKey] || {});
  setPath(nextOverrides, path, parsed);

  // Clean up overrides that match defaults
  const baseTarget = getPath(baseSection, path);
  const overrideTarget = getPath(nextOverrides, path);
  if (overrideTarget === null || deepEqual(overrideTarget, baseTarget)) {
    deletePath(nextOverrides, path);
  }

  if (Object.keys(nextOverrides).length === 0) {
    delete sectionOverrides[sectionKey];
  } else {
    sectionOverrides[sectionKey] = nextOverrides;
  }
  renderConfigCards(defaultSections);
}

function setPath(obj, path, value) {
  if (!path.length) return;
  let cursor = obj;
  for (let i = 0; i < path.length - 1; i += 1) {
    const key = path[i];
    if (typeof cursor[key] !== "object" || cursor[key] === null) {
      cursor[key] = {};
    }
    cursor = cursor[key];
  }
  cursor[path[path.length - 1]] = value;
}

function getPath(obj, path) {
  let cursor = obj;
  for (let i = 0; i < path.length; i += 1) {
    if (cursor == null) return undefined;
    cursor = cursor[path[i]];
  }
  return cursor;
}

function deletePath(obj, path) {
  if (!path.length) return;
  let cursor = obj;
  const stack = [];
  for (let i = 0; i < path.length - 1; i += 1) {
    if (cursor == null) return;
    stack.push([cursor, path[i]]);
    cursor = cursor[path[i]];
  }
  if (cursor && Object.prototype.hasOwnProperty.call(cursor, path[path.length - 1])) {
    delete cursor[path[path.length - 1]];
  }
  // Cleanup empty parents
  for (let i = stack.length - 1; i >= 0; i -= 1) {
    const [parent, key] = stack[i];
    if (parent[key] && typeof parent[key] === "object" && Object.keys(parent[key]).length === 0) {
      delete parent[key];
    }
  }
}

function resetOverrides() {
  sectionOverrides = {};
  renderConfigCards(defaultSections);
  renderPre("config-view", "");
}

function deepMerge(base, overrides) {
  if (!overrides) return cloneDeep(base);
  const result = cloneDeep(base);
  Object.entries(overrides || {}).forEach(([key, value]) => {
    if (value && typeof value === "object" && !Array.isArray(value) && result[key] && typeof result[key] === "object") {
      result[key] = deepMerge(result[key], value);
    } else {
      result[key] = value;
    }
  });
  return result;
}

function buildCurrentConfig() {
  return deepMerge(defaultSections || {}, sectionOverrides || {});
}

async function addRunById(runId) {
  try {
    const status = await fetchJSON(`/runs/${runId}`);
    addRunToBoard(status);
    attachRunInput.value = "";
  } catch (err) {
    alert(err.message);
  }
}

function addRunToBoard(status, meta = {}) {
  renderRunCard(status, meta);
  startPolling(status.run_id);
  switchView("runs-view");
}

function renderRunCard(status, meta = {}) {
  let card = runCards.get(status.run_id);
  if (!card) {
    const panel = document.createElement("section");
    panel.className = "panel";
    panel.dataset.runId = status.run_id;
    panel.innerHTML = `
      <div class="run-card-header">
        <div>
          <p class="eyebrow">${status.workflow.toUpperCase()} • ${status.runtime.toUpperCase()}</p>
          <h3><span class="run-id">${status.run_id}</span> <span class="status-pill status-${status.status}">${status.status}</span></h3>
          <p class="muted run-message"></p>
        </div>
        <div class="run-actions run-card-actions">
          <button class="ghost" data-action="refresh">Refresh</button>
          <button class="ghost" data-action="download-log">Download log</button>
          <button data-action="download-artifacts">Download outputs</button>
          <button class="warn" data-action="cancel">Cancel</button>
          <button class="ghost" data-action="delete">Remove run</button>
        </div>
      </div>
      <div class="run-meta"></div>
      <div class="console log-console">Waiting for log...</div>
    `;
    panel.querySelectorAll("[data-action]").forEach((btn) => {
      btn.addEventListener("click", () => handleRunAction(status.run_id, btn.dataset.action));
    });
    runsContainer.prepend(panel);
    card = {
      panel,
      statusPill: panel.querySelector(".status-pill"),
      messageEl: panel.querySelector(".run-message"),
      metaEl: panel.querySelector(".run-meta"),
      logEl: panel.querySelector(".log-console"),
      timer: null,
      reportType: meta.reportType || null,
      downloaded: false,
    };
    runCards.set(status.run_id, card);
    updateEmptyState();
  }
  updateRunCard(status);
  return card;
}

function updateRunCard(status) {
  const card = runCards.get(status.run_id);
  if (!card) return;
  card.statusPill.textContent = status.status;
  card.statusPill.className = `status-pill status-${status.status}`;
  card.messageEl.textContent = status.message || "";
  if (["completed", "failed", "cancelled"].includes(status.status)) {
    if (card.timer) {
      clearInterval(card.timer);
      card.timer = null;
      updatePollingCount();
    }
    if (!card.downloaded && status.status === "completed") {
      const shouldDownloadZip = card.reportType === 2;
      if (shouldDownloadZip) {
        packageAndDownload(status.run_id);
        card.downloaded = true;
      }
    }
  }
  card.metaEl.innerHTML = [
    ["Status", status.status],
    ["Run dir", status.run_dir],
    ["Log", status.log_path],
    ["PID", status.pid],
    ["Return code", status.return_code],
    ["Scheduler ID", status.scheduler_job_id],
    ["Submitted", formatTime(status.submitted_at)],
    ["Started", formatTime(status.started_at)],
    ["Ended", formatTime(status.ended_at)],
    ["Command", status.command],
  ]
    .filter(([, value]) => value !== null && value !== undefined && value !== "")
    .map(([label, value]) => `<strong>${label}:</strong> ${value}`)
    .join("<br>");
}

async function handleRunAction(runId, action) {
  if (action === "refresh") {
    await refreshRun(runId);
    return;
  }
  if (action === "download-log") {
    await downloadLog(runId);
    return;
  }
  if (action === "download-artifacts") {
    await packageAndDownload(runId);
    return;
  }
  if (action === "cancel") {
    await cancelRun(runId);
    return;
  }
  if (action === "delete") {
    await deleteRun(runId);
  }
}

function startPolling(runId) {
  const card = runCards.get(runId);
  if (!card) return;
  if (card.timer) clearInterval(card.timer);
  const tick = () => refreshRun(runId, { silent: true });
  tick();
  card.timer = setInterval(tick, REFRESH_MS);
  updatePollingCount();
}

function updatePollingCount() {
  const active = Array.from(runCards.values()).filter((card) => card.timer).length;
  pollingCount.textContent = `${active} runs auto-refreshing`;
}

function updateEmptyState() {
  if (!emptyRuns) return;
  emptyRuns.style.display = runCards.size ? "none" : "block";
}

async function refreshRun(runId, opts = {}) {
  try {
    const status = await fetchJSON(`/runs/${runId}`);
    renderRunCard(status);
    await updateLogs(runId);
  } catch (err) {
    if (!opts.silent) alert(err.message);
  }
}

async function updateLogs(runId) {
  const card = runCards.get(runId);
  if (!card) return;
  try {
    const res = await fetchJSON(`/runs/${runId}/logs`);
    const content = res.tail.join("\n");
    card.logEl.textContent = content || "No log output yet.";
    card.logEl.scrollTop = card.logEl.scrollHeight;
  } catch (err) {
    card.logEl.textContent = `Unable to read logs: ${err.message}`;
  }
}

async function cancelRun(runId) {
  try {
    const res = await fetchJSON(`/runs/${runId}/cancel`, { method: "POST" });
    updateRunCard(res);
  } catch (err) {
    alert(err.message);
  }
}

async function packageAndDownload(runId) {
  try {
    await fetchJSON(`/runs/${runId}/artifacts`);
    await downloadFile(`/runs/${runId}/artifacts/download`, `${runId}.tar.gz`);
  } catch (err) {
    alert(err.message);
  }
}

async function downloadLog(runId) {
  try {
    await downloadFile(`/runs/${runId}/logs/download`, `${runId}-snakemake.log`);
  } catch (err) {
    alert(err.message);
  }
}

async function deleteRun(runId) {
  if (!confirm(`Remove run ${runId}? This deletes its run directory and artifacts.`)) return;
  try {
    await fetchJSON(`/runs/${runId}`, { method: "DELETE" });
    const card = runCards.get(runId);
    if (card && card.panel) {
      card.panel.remove();
    }
    if (card && card.timer) {
      clearInterval(card.timer);
    }
    runCards.delete(runId);
    updatePollingCount();
    updateEmptyState();
  } catch (err) {
    alert(err.message);
  }
}

async function uploadAsset(kind) {
  const input = document.getElementById(`${kind === "classifiers" ? "classifier" : "adapter"}-file`);
  if (!input.files.length) {
    alert("Choose a file first");
    return;
  }
  const form = new FormData();
  form.append("file", input.files[0]);
  const res = await fetchJSON(`/${kind}`, { method: "POST", body: form, isForm: true });
  renderPre("run-response", res);
  listAssets();
}

async function listAssets() {
  const classifiers = await fetchJSON("/classifiers");
  renderList("classifier-list", classifiers.items, "classifiers");
  const adapters = await fetchJSON("/adapters");
  renderList("adapter-list", adapters.items, "adapters");
}

function renderList(id, items, kind) {
  const el = document.getElementById(id);
  el.innerHTML = "";
  items.forEach((item) => {
    const li = document.createElement("li");
    li.textContent = item;
    const btn = document.createElement("button");
    btn.textContent = "Delete";
    btn.className = "ghost";
    btn.onclick = async () => {
      await fetchJSON(`/delete/${kind}/${item}`, { method: "POST" });
      listAssets();
    };
    li.appendChild(btn);
    el.appendChild(li);
  });
}

async function fetchJSON(path, opts = {}) {
  const headers = opts.isForm ? {} : { "Content-Type": "application/json" };
  const res = await fetch(`${apiBase}${path}`, { ...opts, headers });
  if (!res.ok) {
    const text = await res.text();
    throw new Error(text || res.statusText);
  }
  return res.json();
}

function renderPre(id, content) {
  const el = document.getElementById(id);
  if (typeof content === "string") {
    el.textContent = content;
  } else {
    el.textContent = JSON.stringify(content, null, 2);
  }
}

async function downloadFile(path, filename) {
  const res = await fetch(`${apiBase}${path}`);
  if (!res.ok) {
    throw new Error(await res.text());
  }
  const blob = await res.blob();
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function formatTime(value) {
  if (!value) return null;
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) return value;
  return date.toLocaleString();
}

async function applyDetectedPaths() {
  try {
    const res = await fetchJSON("/settings/paths");
    const hostRepo = (res.repo_root || "").replace(/\\/g, "/");
    if (hostRepo) {
      autoBindValue = `${hostRepo}:/workspace`;
    }
    if (!cacheDirInput.value) {
      cacheDirInput.value = "/workspace/runtime/cache";
    }
    if (!inputDirInput.value) {
      inputDirInput.value = "/workspace/testing/testing_data";
    }
  } catch (err) {
    console.warn("Could not detect repo path", err);
  }
}

switchView("setup-view");
applyDetectedPaths();
loadDefaultConfig(workflowSelect.value);
listAssets();
