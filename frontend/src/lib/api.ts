import type {
  ArtifactResponse,
  ConfigSchemaResponse,
  ConfigSectionsResponse,
  DeleteRunResponse,
  ListAssetsResponse,
  LogResponse,
  PathsResponse,
  ProfileListResponse,
  RenderConfigResponse,
  RenderConfigWithProfileRequest,
  RunListResponse,
  RunStatus,
  RunSubmissionRequest,
  UploadResponse,
} from "../generated/api";
import {
  assetsListSchema,
  configSchemaResponseSchema,
  configSectionsResponseSchema,
  pathsResponseSchema,
  profileListSchema,
  renderConfigResponseSchema,
  runListSchema,
  runStatusSchema,
  uploadResponseSchema,
} from "./validation";

export class ApiError extends Error {
  status: number;

  constructor(message: string, status: number) {
    super(message);
    this.status = status;
  }
}

const API_BASE = (import.meta.env.VITE_API_BASE || "/api").replace(/\/$/, "");

export function buildApiUrl(path: string): string {
  return `${API_BASE}${path.startsWith("/") ? path : `/${path}`}`;
}

async function requestJson<T>(
  path: string,
  init: RequestInit,
  parser?: { parse: (data: unknown) => T },
): Promise<T> {
  const response = await fetch(buildApiUrl(path), init);
  const contentType = response.headers.get("content-type") || "";
  const isJson = contentType.includes("application/json");
  const payload = isJson ? await response.json().catch(() => null) : await response.text();

  if (!response.ok) {
    const detail =
      payload && typeof payload === "object" && "detail" in payload
        ? String(payload.detail)
        : typeof payload === "string"
          ? payload
          : `Request failed with status ${response.status}`;
    throw new ApiError(detail, response.status);
  }

  return parser ? parser.parse(payload) : (payload as T);
}

async function requestText(path: string): Promise<string> {
  const response = await fetch(buildApiUrl(path), {
    headers: { Accept: "text/plain" },
  });
  const text = await response.text();
  if (!response.ok) {
    throw new ApiError(text || `Request failed with status ${response.status}`, response.status);
  }
  return text;
}

export function getHealth(): Promise<string> {
  return requestText("/health");
}

export function getSettingsPaths(): Promise<PathsResponse> {
  return requestJson("/settings/paths", { method: "GET" }, pathsResponseSchema);
}

export function listRuns(): Promise<RunListResponse> {
  return requestJson("/runs", { method: "GET" }, runListSchema);
}

export function getRun(runId: string): Promise<RunStatus> {
  return requestJson(`/runs/${encodeURIComponent(runId)}`, { method: "GET" }, runStatusSchema);
}

export function submitRun(payload: RunSubmissionRequest): Promise<RunStatus> {
  return requestJson(
    "/runs",
    {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    },
    runStatusSchema,
  );
}

export function getRunLogs(runId: string, lines = 200): Promise<LogResponse> {
  return requestJson(`/runs/${encodeURIComponent(runId)}/logs?lines=${lines}`, { method: "GET" });
}

export function cancelRun(runId: string): Promise<RunStatus> {
  return requestJson(
    `/runs/${encodeURIComponent(runId)}/cancel`,
    { method: "POST" },
    runStatusSchema,
  );
}

export function deleteRun(runId: string): Promise<DeleteRunResponse> {
  return requestJson(`/runs/${encodeURIComponent(runId)}`, { method: "DELETE" });
}

export function packageArtifacts(runId: string): Promise<ArtifactResponse> {
  return requestJson(`/runs/${encodeURIComponent(runId)}/artifacts`, { method: "GET" });
}

export function listProfiles(): Promise<ProfileListResponse> {
  return requestJson("/profiles", { method: "GET" }, profileListSchema);
}

export function getProfile(profileName: string): Promise<Record<string, unknown>> {
  return requestJson(`/profiles/${encodeURIComponent(profileName)}`, { method: "GET" });
}

export function getDefaultConfigSections(workflow: "esv" | "otu"): Promise<ConfigSectionsResponse> {
  return requestJson(
    `/configs/defaults/${encodeURIComponent(workflow)}/sections`,
    { method: "GET" },
    configSectionsResponseSchema,
  );
}

export function renderConfig(
  payload: RenderConfigWithProfileRequest,
): Promise<RenderConfigResponse> {
  return requestJson(
    "/configs/render",
    {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    },
    renderConfigResponseSchema,
  );
}

async function uploadFile(path: string, file: File): Promise<UploadResponse> {
  const formData = new FormData();
  formData.append("file", file);
  return requestJson(path, { method: "POST", body: formData }, uploadResponseSchema);
}

export function listAdapters(): Promise<ListAssetsResponse> {
  return requestJson("/adapters", { method: "GET" }, assetsListSchema);
}

export function uploadAdapter(file: File): Promise<UploadResponse> {
  return uploadFile("/adapters", file);
}

export function deleteAdapter(name: string): Promise<ListAssetsResponse> {
  return requestJson(
    `/adapters/${encodeURIComponent(name)}`,
    { method: "DELETE" },
    assetsListSchema,
  );
}

export function listClassifiers(): Promise<ListAssetsResponse> {
  return requestJson("/classifiers", { method: "GET" }, assetsListSchema);
}

export function uploadClassifier(file: File): Promise<UploadResponse> {
  return uploadFile("/classifiers", file);
}

export function deleteClassifier(name: string): Promise<ListAssetsResponse> {
  return requestJson(
    `/classifiers/${encodeURIComponent(name)}`,
    { method: "DELETE" },
    assetsListSchema,
  );
}

export function getConfigSchema(
  profile: string,
  workflow: string = "esv",
): Promise<ConfigSchemaResponse> {
  return requestJson(
    `/config/schema?profile=${encodeURIComponent(profile)}&workflow=${encodeURIComponent(workflow)}`,
    { method: "GET" },
    configSchemaResponseSchema,
  );
}

