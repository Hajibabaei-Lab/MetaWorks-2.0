/* Generated from the MetaWorks FastAPI schema. Regenerate with `npm run generate:types`. */

export type WorkflowType = "esv" | "otu";
export type RuntimeType = "conda" | "docker" | "apptainer";

export interface ProfileInfo {
  name: string;
  description?: string;
  marker?: string;
  file: string;
}

export interface ProfileListResponse {
  profiles: ProfileInfo[];
}

export interface RenderConfigWithProfileRequest {
  profile?: string;
  workflow?: WorkflowType;
  overrides?: Record<string, unknown>;
}

export interface RunSubmissionRequest {
  profile?: string;
  workflow: WorkflowType;
  run_name: string;
  runtime?: RuntimeType;
  container_image?: string | null;
  bind_paths?: string[];
  cache_dir?: string | null;
  config_overrides?: Record<string, unknown>;
  input_dir: string;
  samples_csv?: string | null;
  sample_source?: "folder" | "csv";
  notes?: string | null;
  cores?: number | null;
  mem_gb?: number | null;
  dry_run?: boolean;
  keep_outputs?: boolean;
}

export interface RunStatus {
  run_id: string;
  workflow: WorkflowType;
  runtime: RuntimeType;
  status: string;
  run_name?: string | null;
  sample_source?: string | null;
  input_dir?: string | null;
  samples_csv?: string | null;
  notes?: string | null;
  config_overrides?: Record<string, unknown> | null;
  scheduler_job_id?: string | null;
  submitted_at?: string | null;
  started_at?: string | null;
  ended_at?: string | null;
  pid?: number | null;
  return_code?: number | null;
  message?: string | null;
  run_dir?: string | null;
  log_path?: string | null;
  command?: string | null;
  artifact_path?: string | null;
  keep_outputs?: boolean | null;
}

export interface RunListResponse {
  runs: RunStatus[];
}

export interface LogResponse {
  run_id: string;
  tail: string[];
}

export interface RenderConfigResponse {
  workflow: WorkflowType;
  merged: string;
}

export interface ConfigSectionsResponse {
  workflow: WorkflowType;
  sections: Record<string, unknown>;
}

export interface UploadResponse {
  name: string;
  path: string;
}

export interface ListAssetsResponse {
  items: string[];
}

export interface ArtifactResponse {
  run_id: string;
  archive_path: string;
}

export interface DeleteRunResponse {
  run_id: string;
  removed_paths: string[];
}

export interface PathsResponse {
  repo_root: string;
  runtime_cache: string;
  allowed_runtimes: string[];
  retention_policy: string;
  default_runtime: string;
  container_image: string;
}

export interface FieldSchema {
  type: string;
  default?: unknown;
  description?: string;
  constraints?: Record<string, unknown>;
  options?: unknown[];
  visible_when?: string;
  label?: string;
  fields?: Record<string, FieldSchema>;
}

export interface SectionSchema {
  label: string;
  enabled_by?: string;
  collapsed?: boolean;
  fields: Record<string, FieldSchema>;
}

export interface ConfigSchemaResponse {
  profile: string;
  marker: string;
  workflow: string;
  sections: Record<string, SectionSchema>;
}
